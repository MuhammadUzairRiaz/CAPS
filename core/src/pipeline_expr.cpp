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

struct Node {
  enum Kind { Num, Var, Neg, Not, Bin, Call } kind = Num;
  double value = 0;
  const std::vector<double>* var = nullptr;
  std::string op;
  std::vector<std::unique_ptr<Node>> args;

  double eval(size_t i) const {
    switch (kind) {
      case Num: return value;
      case Var: return (*var)[i];
      case Neg: return -args[0]->eval(i);
      case Not: return args[0]->eval(i) == 0 ? 1 : 0;
      case Bin: {
        if (op == "&&") return args[0]->eval(i) != 0 && args[1]->eval(i) != 0 ? 1 : 0;
        if (op == "||") return args[0]->eval(i) != 0 || args[1]->eval(i) != 0 ? 1 : 0;
        const double a = args[0]->eval(i), b = args[1]->eval(i);
        if (op == "+") return a + b;
        if (op == "-") return a - b;
        if (op == "*") return a * b;
        if (op == "/") return b != 0 ? a / b : std::nan("");
        if (op == "%") return b != 0 ? std::fmod(a, b) : std::nan("");
        if (op == "^") return std::pow(a, b);
        constexpr double eps = 1e-9;
        if (op == "==") return std::fabs(a - b) <= eps * std::max(1.0, std::fabs(a)) ? 1 : 0;
        if (op == "!=") return std::fabs(a - b) > eps * std::max(1.0, std::fabs(a)) ? 1 : 0;
        if (op == "<") return a < b ? 1 : 0;
        if (op == "<=") return a <= b ? 1 : 0;
        if (op == ">") return a > b ? 1 : 0;
        if (op == ">=") return a >= b ? 1 : 0;
        return 0;
      }
      case Call: {
        const double a = args[0]->eval(i);
        if (op == "abs") return std::fabs(a);
        if (op == "sqrt") return std::sqrt(a);
        if (op == "exp") return std::exp(a);
        if (op == "log") return std::log(a);
        if (op == "floor") return std::floor(a);
        if (op == "ceil") return std::ceil(a);
        if (op == "round") return std::round(a);
        if (op == "min") return std::min(a, args[1]->eval(i));
        if (op == "max") return std::max(a, args[1]->eval(i));
        return 0;
      }
    }
    return 0;
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
      if (p_ < s_.size() && s_[p_] == '(') {
        static const char* one[] = {"abs", "sqrt", "exp", "log", "floor", "ceil", "round"};
        const bool two = name == "min" || name == "max";
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
  if (text.find_first_not_of(" \t\r\n") == std::string::npos) return std::vector<double>(n, 1.0);
  Parser parser(text, lookup);
  const auto root = parser.parse();
  std::vector<double> out(n);
  for (size_t i = 0; i < n; ++i) out[i] = root->eval(i);
  return out;
}

}  // namespace caps
