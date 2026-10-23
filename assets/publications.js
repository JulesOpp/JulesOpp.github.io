// Sort the publication list by year (default), first author or corresponding author.
// Without JavaScript the list stays grouped by year.
(function () {
  var byYear = document.getElementById("pubs-by-year");
  var sorter = document.getElementById("pub-sort");
  if (!byYear || !sorter) return;

  var items = Array.prototype.slice.call(byYear.querySelectorAll("li[data-year]"));
  var views = { year: byYear };

  function list(title, entries) {
    var frag = document.createDocumentFragment();
    var h = document.createElement("h2");
    h.textContent = title + " (" + entries.length + ")";
    var ol = document.createElement("ol");
    ol.className = "pubs";
    ol.reversed = true;
    entries.forEach(function (li) { ol.appendChild(li.cloneNode(true)); });
    frag.appendChild(h);
    frag.appendChild(ol);
    return frag;
  }

  function build(attr, title) {
    var div = document.createElement("div");
    div.hidden = true;
    var yes = items.filter(function (li) { return li.hasAttribute(attr); });
    var no = items.filter(function (li) { return !li.hasAttribute(attr); });
    div.appendChild(list(title, yes));
    div.appendChild(list("Other publications", no));
    byYear.parentNode.insertBefore(div, byYear.nextSibling);
    return div;
  }

  sorter.hidden = false;
  sorter.addEventListener("click", function (event) {
    var button = event.target.closest("button[data-view]");
    if (!button) return;
    var view = button.getAttribute("data-view");
    if (!views[view]) {
      views[view] = view === "first"
        ? build("data-first", "First author")
        : build("data-corr", "Corresponding author");
    }
    Object.keys(views).forEach(function (k) { views[k].hidden = k !== view; });
    sorter.querySelectorAll("button").forEach(function (b) {
      b.setAttribute("aria-pressed", String(b === button));
    });
  });
})();
