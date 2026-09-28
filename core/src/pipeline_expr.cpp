// CAPS pipeline expressions (see caps/pipeline.hpp): a small recursive-descent parser over particle properties,
// evaluated for every particle at once.
#include <cctype>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/pipeline.hpp"

namespace caps {

namespace {

// A value: a number, or a 3-vector (Position, MoleculeCOM, any property with .X .Y .Z components).
struct Val {
  double x = 0, y = 0, z = 0;
  bool vec = false;
  static Val num(double v) { return {v, 0, 0, false}; }
  double scalar() const {
    if (vec) throw std::invalid_argument("a vector where a number is needed: use norm(…) or a component (.X)");
    return x;
  }
};

struct Node {
  enum Kind { Num, Var, VecVar, Neg, Not, Bin, Call } kind = Num;
  double value = 0;
  const std::vector<double>* var = nullptr;
  const std::vector<double>* comp[3] = {nullptr, nullptr, nullptr};
  std::string op;
  std::vector<std::unique_ptr<Node>> args;

  Val eval(size_t i) const {
    switch (kind) {
      case Num: return Val::num(value);
      case Var: return Val::num((*var)[i]);
      case VecVar: return {(*comp[0])[i], (*comp[1])[i], (*comp[2])[i], true};
      case Neg: { Val a = args[0]->eval(i); return {-a.x, -a.y, -a.z, a.vec}; }
      case Not: return Val::num(args[0]->eval(i).scalar() == 0 ? 1 : 0);
      case Bin: {
        if (op == "&&") return Val::num(args[0]->eval(i).scalar() != 0 && args[1]->eval(i).scalar() != 0 ? 1 : 0);
        if (op == "||") return Val::num(args[0]->eval(i).scalar() != 0 || args[1]->eval(i).scalar() != 0 ? 1 : 0);
        const Val A = args[0]->eval(i), B = args[1]->eval(i);
        if (A.vec || B.vec) {   // vectors: + − between vectors, × and ÷ by a number
          if ((op == "+" || op == "-") && A.vec && B.vec) { const double sg = op == "+" ? 1 : -1; return {A.x + sg * B.x, A.y + sg * B.y, A.z + sg * B.z, true}; }
          if (op == "*" && A.vec != B.vec) { const Val& v = A.vec ? A : B; const double k = A.vec ? B.x : A.x; return {v.x * k, v.y * k, v.z * k, true}; }
          if (op == "/" && A.vec && !B.vec) return {A.x / B.x, A.y / B.x, A.z / B.x, true};
          throw std::invalid_argument("'" + op + "' does not apply to these vectors");
        }
        const double a = A.x, b = B.x;
        if (op == "+") return Val::num(a + b);
        if (op == "-") return Val::num(a - b);
        if (op == "*") return Val::num(a * b);
        if (op == "/") return Val::num(b != 0 ? a / b : std::nan(""));
        if (op == "%") return Val::num(b != 0 ? std::fmod(a, b) : std::nan(""));
        if (op == "^") return Val::num(std::pow(a, b));
        constexpr double eps = 1e-9;
        if (op == "==") return Val::num(std::fabs(a - b) <= eps * std::max(1.0, std::fabs(a)) ? 1 : 0);
        if (op == "!=") return Val::num(std::fabs(a - b) > eps * std::max(1.0, std::fabs(a)) ? 1 : 0);
        if (op == "<") return Val::num(a < b ? 1 : 0);
        if (op == "<=") return Val::num(a <= b ? 1 : 0);
        if (op == ">") return Val::num(a > b ? 1 : 0);
        if (op == ">=") return Val::num(a >= b ? 1 : 0);
        return Val::num(0);
      }
      case Call: {
        if (op == "norm") { const Val v = args[0]->eval(i); return Val::num(v.vec ? std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z) : std::fabs(v.x)); }
        if (op == "dot") {
          const Val a = args[0]->eval(i), b = args[1]->eval(i);
          if (!a.vec || !b.vec) throw std::invalid_argument("dot(…) takes two vectors");
          return Val::num(a.x * b.x + a.y * b.y + a.z * b.z);
        }
        const double a = args[0]->eval(i).scalar();
        if (op == "abs") return Val::num(std::fabs(a));
        if (op == "sqrt") return Val::num(std::sqrt(a));
        if (op == "exp") return Val::num(std::exp(a));
        if (op == "log") return Val::num(std::log(a));
        if (op == "floor") return Val::num(std::floor(a));
        if (op == "ceil") return Val::num(std::ceil(a));
        if (op == "round") return Val::num(std::round(a));
        if (op == "min") return Val::num(std::min(a, args[1]->eval(i).scalar()));
        if (op == "max") return Val::num(std::max(a, args[1]->eval(i).scalar()));
        return Val::num(0);
      }
    }
    return Val::num(0);
  }
};

class Parser {
 public:
  Parser(const std::string& text, std::function<const std::vector<double>*(const std::string&)> lookup)
      : s_(text), lookup_(std::move(lookup)) {}

  std::unique_ptr<Node> parse() {
    auto n = expr();
    skip();
    if (p_ < s_.size()) fail("unexpected '" + s_.substr(p_, 1) + "'");
    return n;
  }

 private:
  const std::string& s_;
  size_t p_ = 0;
  std::function<const std::vector<double>*(const std::string&)> lookup_;

  // a name with .X .Y .Z components is a vector (Position, MoleculeCOM …)
  bool vector_name(const std::string& name) { return lookup_(name + ".X") && lookup_(name + ".Y") && lookup_(name + ".Z"); }
  std::unique_ptr<Node> vector_var(const std::string& name, size_t at) {
    auto n = std::make_unique<Node>();
    n->kind = Node::VecVar;
    const char* c[] = {".X", ".Y", ".Z"};
    for (int k = 0; k < 3; ++k)
      if (!(n->comp[k] = lookup_(name + c[k]))) { p_ = at; fail("unknown property " + name); }
    return n;
  }

  [[noreturn]] void fail(const std::string& what) const {
    throw std::invalid_argument(what + " at character " + std::to_string(p_ + 1));
  }
  void skip() { while (p_ < s_.size() && std::isspace(static_cast<unsigned char>(s_[p_]))) ++p_; }
  bool eat(const std::string& t) {
    skip();
    if (s_.compare(p_, t.size(), t) != 0) return false;
    // a word operator must not run into a name (and, or, not)
    if (std::isalpha(static_cast<unsigned char>(t[0])) && p_ + t.size() < s_.size() &&
        (std::isalnum(static_cast<unsigned char>(s_[p_ + t.size()])) || s_[p_ + t.size()] == '_'))
      return false;
    p_ += t.size();
    return true;
  }
  static std::unique_ptr<Node> bin(std::string op, std::unique_ptr<Node> a, std::unique_ptr<Node> b) {
    auto n = std::make_unique<Node>();
    n->kind = Node::Bin;
    n->op = std::move(op);
    n->args.push_back(std::move(a));
    n->args.push_back(std::move(b));
    return n;
  }

  std::unique_ptr<Node> expr() { return orr(); }
  std::unique_ptr<Node> orr() {
    auto a = andd();
    while (eat("||") || eat("or")) a = bin("||", std::move(a), andd());
    return a;
  }
  std::unique_ptr<Node> andd() {
    auto a = cmp();
    while (eat("&&") || eat("and")) a = bin("&&", std::move(a), cmp());
    return a;
  }
  std::unique_ptr<Node> cmp() {
    auto a = add();
    for (const char* op : {"==", "!=", "<=", ">=", "<", ">"})
      if (eat(op)) return bin(op, std::move(a), add());
    if (eat("=")) return bin("==", std::move(a), add());
    return a;
  }
  std::unique_ptr<Node> add() {
    auto a = mul();
    for (;;) {
      if (eat("+")) a = bin("+", std::move(a), mul());
      else if (eat("-")) a = bin("-", std::move(a), mul());
      else return a;
    }
  }
  std::unique_ptr<Node> mul() {
    auto a = unary();
    for (;;) {
      if (eat("*")) a = bin("*", std::move(a), unary());
      else if (eat("/")) a = bin("/", std::move(a), unary());
      else if (eat("%")) a = bin("%", std::move(a), unary());
      else return a;
    }
  }
  std::unique_ptr<Node> unary() {
    if (eat("-")) {
      auto n = std::make_unique<Node>();
      n->kind = Node::Neg;
      n->args.push_back(unary());
      return n;
    }
    if (eat("+")) return unary();
    skip();
    if ((p_ < s_.size() && s_[p_] == '!' && (p_ + 1 >= s_.size() || s_[p_ + 1] != '=') && (++p_, true)) || eat("not")) {
      auto n = std::make_unique<Node>();
      n->kind = Node::Not;
      n->args.push_back(unary());
      return n;
    }
    auto a = primary();
    if (eat("^")) return bin("^", std::move(a), unary());
    return a;
  }
  std::unique_ptr<Node> primary() {
    skip();
    if (p_ >= s_.size()) fail("expression ends early");
    const char c = s_[p_];
    if (c == '(') {
      ++p_;
      auto n = expr();
      if (!eat(")")) fail("missing ')'");
      return n;
    }
    if (c == '"' || c == '\'') {
      const size_t e = s_.find(c, p_ + 1);
      if (e == std::string::npos) fail("unterminated string");
      const std::string sym = s_.substr(p_ + 1, e - p_ - 1);
      const int z = element_from_symbol(sym);
      if (z == 0) fail("\"" + sym + "\" is not an element");
      p_ = e + 1;
      auto n = std::make_unique<Node>();
      n->value = z;
      return n;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      size_t used = 0;
      double v = 0;
      try { v = std::stod(s_.substr(p_), &used); } catch (...) { fail("bad number"); }
      p_ += used;
      auto n = std::make_unique<Node>();
      n->value = v;
      return n;
    }
    if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
      const size_t b = p_;
      while (p_ < s_.size() && (std::isalnum(static_cast<unsigned char>(s_[p_])) || s_[p_] == '_' || s_[p_] == '.')) ++p_;
      const std::string name = s_.substr(b, p_ - b);
      skip();
      if (p_ < s_.size() && s_[p_] == '(' && vector_name(name)) {   // MoleculeCOM(MoleculeIdentifier): the argument names the grouping
        ++p_;
        int depth = 1;
        while (p_ < s_.size() && depth > 0) depth += s_[p_] == '(' ? 1 : s_[p_] == ')' ? -1 : 0, ++p_;
        if (depth != 0) fail("missing ')'");
        return vector_var(name, b);
      }
      if (p_ < s_.size() && s_[p_] == '(') {
        static const char* one[] = {"abs", "sqrt", "exp", "log", "floor", "ceil", "round", "norm"};
        const bool two = name == "min" || name == "max" || name == "dot";
        bool known = two;
        for (const char* f : one) known |= name == f;
        if (!known) { p_ = b; fail("unknown function " + name); }
        ++p_;
        auto n = std::make_unique<Node>();
        n->kind = Node::Call;
        n->op = name;
        n->args.push_back(expr());
        if (two) {
          if (!eat(",")) fail(name + " takes two values");
          n->args.push_back(expr());
        }
        if (!eat(")")) fail("missing ')'");
        return n;
      }
      if (name == "pi") { auto n = std::make_unique<Node>(); n->value = M_PI; return n; }
      if (name == "true") { auto n = std::make_unique<Node>(); n->value = 1; return n; }
      if (name == "false") { auto n = std::make_unique<Node>(); n->value = 0; return n; }
      const auto* v = lookup_(name);
      if (!v && vector_name(name)) return vector_var(name, b);
      if (!v) { p_ = b; fail("unknown property " + name); }
      auto n = std::make_unique<Node>();
      n->kind = Node::Var;
      n->var = v;
      return n;
    }
    fail("unexpected '" + std::string(1, c) + "'");
  }
};

}  // namespace

std::vector<double> evaluate_expression(const PipelineState& st, const std::string& expr) {
  const size_t n = st.system.atoms.size();
  std::map<std::string, std::vector<double>> cache;
  auto lookup = [&](const std::string& name) -> const std::vector<double>* {
    auto it = cache.find(name);
    if (it != cache.end()) return &it->second;
    std::vector<double> v;
    if (!property_values(st, name, v)) return nullptr;
    return &(cache[name] = std::move(v));
  };
  std::string text = expr;
  for (const auto& [from, to] : {std::pair<std::string, std::string>{"\u2212", "-"}, {"\u00D7", "*"}, {"\u00B7", "*"}}) {   // − × ·
    for (size_t k = text.find(from); k != std::string::npos; k = text.find(from, k + to.size())) text.replace(k, from.size(), to);
  }
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) return std::vector<double>(n, 1.0);
  Parser parser(text, lookup);
  const auto root = parser.parse();
  std::vector<double> out(n);
  for (size_t i = 0; i < n; ++i) out[i] = root->eval(i).scalar();
  return out;
}

std::vector<double> evaluate_pair_expression(const PipelineState& st, const std::string& expr, const std::vector<uint32_t>& neighbour,
                                             const std::vector<Vec3>& delta) {
  const size_t np = neighbour.size();
  std::map<std::string, std::vector<double>> cache;
  auto lookup = [&](const std::string& name) -> const std::vector<double>* {
    auto it = cache.find(name);
    if (it != cache.end()) return &it->second;
    std::vector<double> v(np);
    if (name == "Distance") {
      for (size_t k = 0; k < np; ++k) v[k] = norm(delta[k]);
    } else if (name == "Delta.X" || name == "Delta.Y" || name == "Delta.Z") {
      const size_t c = size_t(name.back() - 'X');
      for (size_t k = 0; k < np; ++k) v[k] = delta[k][c];
    } else {   // the neighbour's own property
      std::vector<double> per;
      if (!property_values(st, name, per)) return nullptr;
      for (size_t k = 0; k < np; ++k) v[k] = per[neighbour[k]];
    }
    return &(cache[name] = std::move(v));
  };
  std::string text = expr;
  for (const auto& [from, to] : {std::pair<std::string, std::string>{"\u2212", "-"}, {"\u00D7", "*"}, {"\u00B7", "*"}}) {
    for (size_t k = text.find(from); k != std::string::npos; k = text.find(from, k + to.size())) text.replace(k, from.size(), to);
  }
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) return std::vector<double>(np, 1.0);
  Parser parser(text, lookup);
  const auto root = parser.parse();
  std::vector<double> out(np);
  for (size_t k = 0; k < np; ++k) out[k] = root->eval(k).scalar();
  return out;
}

}  // namespace caps
