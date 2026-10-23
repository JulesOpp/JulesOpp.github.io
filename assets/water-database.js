// Water Database: builds the tables, the capacity-vs-humidity plot and the CSV downloads from
// data/water-data.js, so editing the data there updates the whole page. Works served or opened from disk.
(function () {
  var SVG = "http://www.w3.org/2000/svg";

  // ---- CSV and table helpers

  function parseCSV(text) {
    var rows = [], row = [], cell = "", quoted = false;
    for (var i = 0; i < text.length; i++) {
      var c = text[i];
      if (quoted) {
        if (c === '"' && text[i + 1] === '"') { cell += '"'; i++; }
        else if (c === '"') quoted = false;
        else cell += c;
      } else if (c === '"') quoted = true;
      else if (c === ",") { row.push(cell); cell = ""; }
      else if (c === "\n" || c === "\r") {
        if (c === "\r" && text[i + 1] === "\n") i++;
        row.push(cell); rows.push(row); row = []; cell = "";
      } else cell += c;
    }
    if (cell || row.length) { row.push(cell); rows.push(row); }
    while (rows.length && rows[rows.length - 1].every(function (x) { return !x.trim(); })) rows.pop();
    return rows;
  }

  function cellNode(text) {
    var t = text.trim(), url = null, label = t;
    if (/^https?:\/\//.test(t)) url = t;
    else if (/^(dx\.)?doi\.org\//.test(t)) url = "https://" + t;
    else if (/^10\.\d{4,}\/\S+$/.test(t)) url = "https://doi.org/" + t;
    if (!url) return document.createTextNode(t);
    var a = document.createElement("a");
    a.href = url;
    a.textContent = label.replace(/^https?:\/\/(dx\.)?(doi\.org\/)?/, "");
    return a;
  }

  function renderTable(rows, container, hasFooter) {
    var table = document.createElement("table");
    var head = table.createTHead().insertRow();
    rows[0].forEach(function (h) {
      var th = document.createElement("th");
      th.textContent = h;
      head.appendChild(th);
    });
    var body = rows.slice(1), foot = null;
    if (hasFooter) foot = body.pop();
    var tbody = table.createTBody();
    body.forEach(function (r) {
      var tr = tbody.insertRow();
      r.forEach(function (c) { tr.insertCell().appendChild(cellNode(c)); });
    });
    if (foot) {
      var tf = table.createTFoot().insertRow();
      foot.forEach(function (c) { tf.insertCell().textContent = c; });
    }
    container.textContent = "";
    container.appendChild(table);
  }

  // ---- Plot

  function el(name, attrs, parent) {
    var e = document.createElementNS(SVG, name);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  // Bounding curve from the Isoreticular Curves paper (Figure 5): Qmax = -0.25 (1/ln a)^2 - 1.75 (1/ln a).
  // Drawn up to Qmax = 3.0 g/g, which it reaches at 1/ln a = -3, i.e. a = exp(-1/3) = 71.7% RH.
  var BOUND_MAX_Q = 3.0, BOUND_END_RH = 100 * Math.exp(-1 / 3);
  function bound(a) { var u = 1 / Math.log(a); return -0.25 * u * u - 1.75 * u; }

  function toPoints(rows, series, nameCol, refFor) {
    var pts = [];
    rows.slice(1).forEach(function (r) {
      var a = parseFloat(r[1]), q = parseFloat(r[2]);
      if (!r[nameCol] || !r[nameCol].trim() || isNaN(a) || isNaN(q)) return;
      var rh = a > 1 ? a : a * 100; // the CSVs store the step as a fraction (P/P0)
      pts.push({ name: r[nameCol].trim(), rh: rh, q: q, series: series(r), ref: refFor(r) });
    });
    return pts;
  }

  function drawPlot(figure, points) {
    var narrow = figure.clientWidth < 520; // phones get a smaller canvas so labels stay legible
    var W = narrow ? 400 : 640, H = narrow ? 330 : 420, m = { l: narrow ? 48 : 56, r: 16, t: 12, b: 46 };
    var yMax = Math.max(BOUND_MAX_Q, Math.ceil(Math.max.apply(null, points.map(function (p) { return p.q; })) * 2) / 2);
    var x = function (v) { return m.l + (v / 100) * (W - m.l - m.r); };
    var y = function (v) { return H - m.b - (v / yMax) * (H - m.t - m.b); };

    var svg = el("svg", { viewBox: "0 0 " + W + " " + H, class: narrow ? "plot narrow" : "plot", role: "img",
      "aria-label": "Maximum gravimetric capacity versus relative humidity of the water uptake step for " + points.length + " MOFs and COFs" });
    el("rect", { x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, class: "plot-frame" }, svg);

    // Region above the bounding curve (high capacity at low humidity), which no sorbent reaches
    var line = [x(0) + "," + y(0)];
    for (var rh = 0.5; rh < BOUND_END_RH; rh += 0.5) line.push(x(rh).toFixed(1) + "," + y(bound(rh / 100)).toFixed(1));
    line.push(x(BOUND_END_RH).toFixed(1) + "," + y(BOUND_MAX_Q).toFixed(1));
    el("path", { d: "M" + line.join(" L") + " L" + x(BOUND_END_RH).toFixed(1) + "," + y(yMax) + " L" + x(0) + "," + y(yMax) + " Z", class: "plot-forbidden" }, svg);
    el("polyline", { points: line.join(" "), class: "plot-bound" }, svg);

    // Axes
    for (var t = 0; t <= 100; t += 20) {
      el("line", { x1: x(t), x2: x(t), y1: y(0), y2: y(0) + 5, class: "plot-tick" }, svg);
      el("text", { x: x(t), y: y(0) + 20, "text-anchor": "middle", class: "plot-label" }, svg).textContent = t;
    }
    for (var v = 0; v <= yMax + 1e-9; v += 0.5) {
      el("line", { x1: m.l - 5, x2: m.l, y1: y(v), y2: y(v), class: "plot-tick" }, svg);
      el("text", { x: m.l - 9, y: y(v) + 4, "text-anchor": "end", class: "plot-label" }, svg).textContent = v.toFixed(1);
    }
    el("text", { x: (m.l + W - m.r) / 2, y: H - 8, "text-anchor": "middle", class: "plot-axis" }, svg).textContent = "Relative humidity of step (%)";
    el("text", { x: -(m.t + H - m.b) / 2, y: 14, transform: "rotate(-90)", "text-anchor": "middle", class: "plot-axis" }, svg).textContent = "Max capacity (g/g)";

    // Points, one group per series
    var shapes = {
      circle: function (g, px, py) { return el("circle", { cx: px, cy: py, r: 3.6 }, g); },
      square: function (g, px, py) { return el("rect", { x: px - 3.6, y: py - 3.6, width: 7.2, height: 7.2 }, g); }
    };
    var groups = {};
    SERIES.forEach(function (s) { groups[s.id] = el("g", { class: "plot-series " + s.id }, svg); });
    points.forEach(function (p) {
      var s = SERIES_BY_ID[p.series];
      var node = shapes[s.shape](groups[s.id], x(p.rh).toFixed(1), y(p.q).toFixed(1));
      el("title", {}, node).textContent = p.name + ": " + fmt(p.rh) + "% RH, " + p.q + " g/g";
      p.px = x(p.rh); p.py = y(p.q);
    });
    var marker = el("circle", { r: 8, class: "plot-marker", visibility: "hidden" }, svg);

    figure.textContent = "";
    figure.appendChild(svg);

    // Legend that toggles series, and a readout for the nearest point (works for touch too)
    var legend = document.createElement("div");
    legend.className = "plot-legend";
    SERIES.forEach(function (s) {
      var n = points.filter(function (p) { return p.series === s.id; }).length;
      var b = document.createElement("button");
      b.type = "button";
      b.className = "legend-item " + s.id;
      b.setAttribute("aria-pressed", "true");
      b.innerHTML = '<svg viewBox="0 0 12 12" aria-hidden="true">' + s.icon + "</svg>";
      b.appendChild(document.createTextNode(s.label + " (" + n + ")"));
      b.addEventListener("click", function () {
        var on = b.getAttribute("aria-pressed") !== "true";
        b.setAttribute("aria-pressed", String(on));
        groups[s.id].style.display = on ? "" : "none";
        if (!on && pinned && pinned.series === s.id) unpin();
      });
      legend.appendChild(b);
    });
    var bl = document.createElement("span");
    bl.className = "legend-item bound";
    bl.innerHTML = '<svg viewBox="0 0 12 12" aria-hidden="true"><path d="M1 11 Q6 9 11 1"/></svg>Bounding curve';
    legend.appendChild(bl);
    figure.appendChild(legend);

    var readout = document.createElement("p");
    readout.className = "plot-readout";
    readout.setAttribute("aria-live", "polite");
    var HINT = "Hover over a point to see the material; click or tap it to keep it selected.";
    readout.textContent = HINT;
    figure.appendChild(readout);

    // Hovering previews a point; clicking pins it until the same point is clicked again.
    var pinned = null;
    function nearest(evt) {
      var r = svg.getBoundingClientRect(), sx = W / r.width;
      var mx = (evt.clientX - r.left) * sx, my = (evt.clientY - r.top) * sx, best = null, bd = 14 * sx;
      points.forEach(function (p) {
        if (groups[p.series].style.display === "none") return;
        var d = Math.hypot(p.px - mx, p.py - my);
        if (d < bd) { bd = d; best = p; }
      });
      return best;
    }
    function show(p) {
      marker.setAttribute("cx", p.px); marker.setAttribute("cy", p.py); marker.setAttribute("visibility", "visible");
      marker.classList.toggle("pinned", p === pinned);
      readout.textContent = "";
      var b = document.createElement("strong"); b.textContent = p.name;
      readout.appendChild(b);
      readout.appendChild(document.createTextNode(": step at " + fmt(p.rh) + "% RH, " + p.q + " g/g · " + SERIES_BY_ID[p.series].label + (p.ref ? " · " : "")));
      if (p.ref) readout.appendChild(cellNode(p.ref));
      if (p === pinned) readout.appendChild(document.createTextNode(" (click again to unselect)"));
    }
    function unpin() {
      pinned = null;
      marker.setAttribute("visibility", "hidden");
      readout.textContent = HINT;
    }
    svg.addEventListener("pointermove", function (evt) {
      if (pinned || evt.pointerType !== "mouse") return;
      var p = nearest(evt);
      if (p) show(p); else { marker.setAttribute("visibility", "hidden"); readout.textContent = HINT; }
    });
    svg.addEventListener("click", function (evt) {
      var p = nearest(evt);
      if (!p) return;
      if (p === pinned) { unpin(); return; }
      pinned = p;
      show(p);
    });
  }

  function fmt(v) { return Math.round(v * 10) / 10; }

  var SERIES = [
    { id: "s-mof", label: "MOFs", shape: "circle", icon: '<circle cx="6" cy="6" r="4"/>' },
    { id: "s-cof", label: "COFs", shape: "square", icon: '<rect x="2" y="2" width="8" height="8"/>' }
  ];
  var SERIES_BY_ID = {};
  SERIES.forEach(function (s) { SERIES_BY_ID[s.id] = s; });

  // ---- Load

  var DATA = window.WATER_DATA;
  if (!DATA) {
    document.querySelectorAll(".needs-data").forEach(function (n) { n.hidden = false; });
    return;
  }
  var csv = {};
  Object.keys(DATA).forEach(function (k) { csv[k] = DATA[k].replace(/^\n/, ""); });

  document.querySelectorAll("a[data-csv]").forEach(function (a) {
    var blob = new Blob(["\ufeff" + csv[a.getAttribute("data-csv")] + "\n"], { type: "text/csv;charset=utf-8" });
    a.href = URL.createObjectURL(blob);
  });

  var mof = parseCSV(csv.mof), cof = parseCSV(csv.cof), reg = parseCSV(csv.regression);
  renderTable(mof, document.getElementById("table-mof"), false);
  renderTable(cof, document.getElementById("table-cof"), false);
  renderTable(reg.slice(1), document.getElementById("table-regression"), true); // row 1 is a fit summary
  var isReview = function (r) { return /^[\d,\s]+$/.test(r[3].trim()); }; // numbered refs point into the Chem. Rev. compilation
  var points = toPoints(mof, function () { return "s-mof"; },
      0, function (r) { return isReview(r) ? "10.1021/acs.chemrev.9b00746" : r[3]; })
    .concat(toPoints(cof, function () { return "s-cof"; }, 0, function (r) { return r[5]; }));
  drawPlot(document.getElementById("water-plot"), points);
  var named = function (rows) { return rows.slice(1).filter(function (r) { return r[0].trim(); }).length; };
  document.querySelectorAll("[data-count]").forEach(function (n) {
    var k = n.getAttribute("data-count");
    n.textContent = k === "mof" ? named(mof) : k === "cof" ? named(cof) : points.length;
  });
})();
