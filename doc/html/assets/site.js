// InstanceSpace documentation: search, expand/collapse, copy buttons,
// contents drawer and "On this page" highlighting. Plain JavaScript with
// no network access, so it also works inside MATLAB's Help browser.
(function () {
  "use strict";

  // Open a <details> that holds the URL's #target, e.g. #arg-opts.
  function revealHash() {
    if (!location.hash) return;
    var el = document.getElementById(decodeURIComponent(location.hash.slice(1)));
    for (var n = el; n; n = n.parentElement) {
      if (n.tagName === "DETAILS") n.open = true;
    }
    if (el) el.scrollIntoView();
  }

  function setupToggles() {
    document.querySelectorAll(".toggle-all").forEach(function (btn) {
      var section = btn.closest("section");
      var kind = btn.getAttribute("data-target");
      var items = section.querySelectorAll(":scope > details." + kind);
      function sync() {
        var allOpen = Array.prototype.every.call(items, function (d) { return d.open; });
        btn.textContent = allOpen ? "Collapse all" : "Expand all";
      }
      btn.addEventListener("click", function () {
        var open = btn.textContent === "Expand all";
        section.querySelectorAll("details." + kind).forEach(function (d) { d.open = open; });
        sync();
      });
      items.forEach(function (d) { d.addEventListener("toggle", sync); });
      sync();
    });
  }

  function setupCopy() {
    document.querySelectorAll(".code .copy").forEach(function (btn) {
      btn.addEventListener("click", function () {
        var text = btn.parentElement.querySelector("pre").innerText;
        var done = function () {
          btn.textContent = "Copied";
          setTimeout(function () { btn.textContent = "Copy"; }, 1500);
        };
        if (navigator.clipboard && window.isSecureContext) {
          navigator.clipboard.writeText(text).then(done, function () {});
        } else {
          var ta = document.createElement("textarea");
          ta.value = text; document.body.appendChild(ta); ta.select();
          try { document.execCommand("copy"); done(); } catch (e) { /* no clipboard */ }
          document.body.removeChild(ta);
        }
      });
    });
  }

  function setupMenu() {
    var btn = document.querySelector(".topbar .menu");
    if (!btn) return;
    btn.addEventListener("click", function () {
      var open = document.body.classList.toggle("toc-open");
      btn.setAttribute("aria-expanded", open ? "true" : "false");
    });
    document.querySelector(".content").addEventListener("click", function () {
      document.body.classList.remove("toc-open");
      btn.setAttribute("aria-expanded", "false");
    });
    var current = document.querySelector(".toc a.current");
    if (current) current.scrollIntoView({ block: "center" });
  }

  function setupOnPage() {
    var links = document.querySelectorAll(".onpage a");
    if (!links.length || !("IntersectionObserver" in window)) return;
    var byId = {};
    links.forEach(function (a) { byId[a.getAttribute("href").slice(1)] = a; });
    var observer = new IntersectionObserver(function (entries) {
      entries.forEach(function (e) {
        if (e.isIntersecting && byId[e.target.id]) {
          links.forEach(function (a) { a.classList.remove("active"); });
          byId[e.target.id].classList.add("active");
        }
      });
    }, { rootMargin: "-60px 0px -70% 0px" });
    document.querySelectorAll("section.refsect").forEach(function (s) { observer.observe(s); });
  }

  // Rank pages: title match > heading match > purpose > body text.
  function search(query) {
    var index = window.ISA_SEARCH || [];
    var terms = query.toLowerCase().split(/\s+/).filter(Boolean);
    var results = [];
    index.forEach(function (page) {
      var score = 0, title = page.t.toLowerCase(), purpose = (page.p || "").toLowerCase();
      var heads = page.h.join(" ").toLowerCase(), body = page.x.toLowerCase();
      var hit = terms.every(function (t) {
        var s = 0;
        if (title === t) s += 100;
        else if (title.indexOf(t) === 0) s += 40;
        else if (title.indexOf(t) >= 0) s += 25;
        if (heads.indexOf(t) >= 0) s += 10;
        if (purpose.indexOf(t) >= 0) s += 6;
        if (body.indexOf(t) >= 0) s += 2;
        score += s;
        return s > 0;
      });
      if (hit) results.push({ page: page, score: score });
    });
    results.sort(function (a, b) { return b.score - a.score || a.page.t.localeCompare(b.page.t); });
    return results.slice(0, 12);
  }

  function setupSearch() {
    var input = document.getElementById("search");
    var list = document.getElementById("search-results");
    if (!input || !list) return;
    var active = -1;
    function render() {
      var q = input.value.trim();
      list.innerHTML = "";
      active = -1;
      if (!q) { list.hidden = true; return; }
      var results = search(q);
      if (!results.length) {
        list.innerHTML = '<li class="r-none">No matches</li>';
      }
      results.forEach(function (r) {
        var li = document.createElement("li");
        var a = document.createElement("a");
        a.href = r.page.u;
        var t = document.createElement("span");
        t.className = "r-title"; t.textContent = r.page.t;
        var s = document.createElement("span");
        s.className = "r-sub"; s.textContent = r.page.p || "";
        a.appendChild(t); a.appendChild(s); li.appendChild(a); list.appendChild(li);
      });
      list.hidden = false;
    }
    function move(step) {
      var items = list.querySelectorAll("a");
      if (!items.length) return;
      if (active >= 0) items[active].classList.remove("active");
      active = (active + step + items.length) % items.length;
      items[active].classList.add("active");
      items[active].scrollIntoView({ block: "nearest" });
    }
    input.addEventListener("input", render);
    input.addEventListener("keydown", function (e) {
      if (e.key === "ArrowDown") { move(1); e.preventDefault(); }
      else if (e.key === "ArrowUp") { move(-1); e.preventDefault(); }
      else if (e.key === "Enter") {
        var items = list.querySelectorAll("a");
        var target = items[active >= 0 ? active : 0];
        if (target) location.href = target.href;
      } else if (e.key === "Escape") { input.value = ""; render(); }
    });
    document.addEventListener("click", function (e) {
      if (!e.target.closest(".search")) list.hidden = true;
    });
    document.addEventListener("keydown", function (e) {
      if (e.key === "/" && document.activeElement !== input) { input.focus(); e.preventDefault(); }
    });
  }

  document.addEventListener("DOMContentLoaded", function () {
    setupToggles();
    setupCopy();
    setupMenu();
    setupOnPage();
    setupSearch();
    revealHash();
  });
  window.addEventListener("hashchange", revealHash);
})();
