// CAPS documentation: copy buttons, shell prompts, maths, the reference filter, the menu and the table of contents
(function () {
  "use strict";
  function ready(f) { if (document.readyState !== "loading") f(); else document.addEventListener("DOMContentLoaded", f); }
  ready(function () {
    // shell blocks: comments dimmed
    document.querySelectorAll('.code[data-lang="shell"] code').forEach(function (c) {
      c.innerHTML = c.innerHTML.split("\n").map(function (l) {
        return /^\s*#/.test(l) ? '<span class="c-com">' + l + "</span>" : l;
      }).join("\n");
    });
    document.querySelectorAll(".code").forEach(function (box) {
      var b = document.createElement("button");
      b.className = "copy"; b.type = "button"; b.textContent = "copy";
      b.addEventListener("click", function () {
        var t = box.querySelector("code").innerText;
        var done = function () { b.textContent = "copied"; setTimeout(function () { b.textContent = "copy"; }, 1400); };
        if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(t).then(done, function () { select(box); });
        else select(box);
      });
      box.appendChild(b);
    });
    function select(box) {
      var r = document.createRange(); r.selectNodeContents(box.querySelector("code"));
      var s = window.getSelection(); s.removeAllRanges(); s.addRange(r);
    }
    // maths (KaTeX loads deferred on the pages that have equations)
    if (window.renderMathInElement) {
      window.renderMathInElement(document.querySelector("main"), {
        delimiters: [{ left: "$$", right: "$$", display: true }, { left: "\\(", right: "\\)", display: false }],
        throwOnError: false
      });
    }
    // the menu on narrow screens
    var menu = document.querySelector(".menu"), links = document.getElementById("nav");
    if (menu && links) menu.addEventListener("click", function () {
      var open = links.classList.toggle("open"); menu.setAttribute("aria-expanded", open ? "true" : "false");
    });
    var toc = document.querySelector("nav.toc");
    if (toc) {
      var h = toc.querySelector("h2");
      if (h) h.addEventListener("click", function () { toc.classList.toggle("open"); });
      // the section in view
      var links2 = Array.prototype.slice.call(toc.querySelectorAll('a[href^="#"]'));
      var targets = links2.map(function (a) { return document.getElementById(decodeURIComponent(a.getAttribute("href").slice(1))); });
      if ("IntersectionObserver" in window) {
        var io = new IntersectionObserver(function (es) {
          es.forEach(function (e) {
            if (!e.isIntersecting) return;
            var i = targets.indexOf(e.target);
            links2.forEach(function (a, k) { a.classList.toggle("on", k === i); });
          });
        }, { rootMargin: "-80px 0px -70% 0px" });
        targets.forEach(function (t) { if (t) io.observe(t); });
      }
    }
    // the command reference's filter
    var f = document.getElementById("cmdfilter");
    if (f) f.addEventListener("input", function () {
      var q = f.value.trim().toLowerCase();
      document.querySelectorAll("article.cmd").forEach(function (a) {
        a.hidden = q && a.textContent.toLowerCase().indexOf(q) < 0;
      });
      document.querySelectorAll("section.group").forEach(function (s) {
        s.hidden = !!q && !s.querySelector("article.cmd:not([hidden])");
      });
    });
  });
})();
// modules: the CLI / Python tabs and full-size screenshots
(function () {
  "use strict";
  function ready(f) { if (document.readyState !== "loading") f(); else document.addEventListener("DOMContentLoaded", f); }
  ready(function () {
    document.querySelectorAll(".tabs").forEach(function (bar) {
      var panels = [], n = bar.nextElementSibling;
      while (n && n.classList.contains("tabpanel")) { panels.push(n); n = n.nextElementSibling; }
      bar.querySelectorAll("button").forEach(function (b, i) {
        b.addEventListener("click", function () {
          bar.querySelectorAll("button").forEach(function (x, k) { x.setAttribute("aria-selected", k === i ? "true" : "false"); });
          panels.forEach(function (p, k) { p.hidden = k !== i; });
        });
      });
    });
    var zooms = document.querySelectorAll("figure.shot .zoom");
    if (!zooms.length || typeof HTMLDialogElement === "undefined") return;
    var dlg = document.createElement("dialog");
    dlg.className = "lightbox";
    dlg.innerHTML = '<img alt="">';
    document.body.appendChild(dlg);
    dlg.addEventListener("click", function () { dlg.close(); });
    zooms.forEach(function (z) {
      z.addEventListener("click", function () {
        var img = z.querySelector("img"), big = dlg.querySelector("img");
        big.src = img.src; big.alt = img.alt;
        dlg.showModal();
      });
    });
  });
})();
