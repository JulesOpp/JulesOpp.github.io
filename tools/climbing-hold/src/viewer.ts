// A small WebGL viewer, adapted from the vial-holder post: one shaded
// mesh, a 10 mm grid on the wall plane (x = 0), orbit controls.
//
// The wall can lean. The angle is measured past vertical, as gyms do:
// 0 = vertical, positive = overhang, negative = slab, 90 = roof. The hold
// frame (x off the wall, z up the wall) is turned into the world by
// W = Ry(angle); the floor stays horizontal in the world.

const VS = `attribute vec3 aP; attribute vec3 aN; uniform mat4 uM; uniform mat4 uW;
varying vec3 vN; varying vec3 vW;
void main() { vec4 w = uW * vec4(aP, 1.0); vN = mat3(uW) * aN; vW = w.xyz; gl_Position = uM * w; }`;
const FS = `#ifdef GL_FRAGMENT_PRECISION_HIGH
precision highp float;
#else
precision mediump float;
#endif
varying vec3 vN; varying vec3 vW;
uniform vec3 uCol; uniform vec3 uEye; uniform vec3 uKey;
void main() {
  vec3 n = normalize(vN), v = normalize(uEye - vW), l = normalize(uKey);
  if (dot(n, v) < 0.0) n = -n;
  float d1 = max(dot(n, l), 0.0), d2 = max(dot(n, v), 0.0);
  float sp = pow(max(dot(n, normalize(l + v)), 0.0), 40.0);
  float rim = pow(1.0 - d2, 3.0);
  vec3 c = uCol * (0.2 + 0.6 * d1 + 0.3 * d2) + vec3(0.14) * sp + vec3(0.1) * rim;
  gl_FragColor = vec4(c, 1.0);
}`;
const LVS = `attribute vec3 aP; uniform mat4 uM; void main() { gl_Position = uM * vec4(aP, 1.0); }`;
const LFS = `precision mediump float; uniform vec4 uCol; void main() { gl_FragColor = uCol; }`;

type Mat = ArrayLike<number>;

function compile(gl: WebGLRenderingContext, vs: string, fs: string) {
  const prog = gl.createProgram()!;
  for (const [type, src] of [[gl.VERTEX_SHADER, vs], [gl.FRAGMENT_SHADER, fs]] as const) {
    const s = gl.createShader(type)!;
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s) ?? "shader");
    gl.attachShader(prog, s);
  }
  gl.linkProgram(prog);
  if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog) ?? "link");
  return prog;
}

function perspective(fy: number, asp: number, near: number, far: number): number[] {
  const f = 1 / Math.tan(fy / 2), nf = 1 / (near - far);
  return [f / asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0];
}
function ortho(hh: number, asp: number, near: number, far: number): number[] {
  const nf = 1 / (near - far);
  return [1 / (hh * asp), 0, 0, 0, 0, 1 / hh, 0, 0, 0, 0, 2 * nf, 0, 0, 0, (far + near) * nf, 1];
}
/** Rotation about +y: hold x (off the wall) -> (cos, 0, -sin). */
function rotY(a: number): number[] {
  const c = Math.cos(a), s = Math.sin(a);
  return [c, 0, -s, 0, 0, 1, 0, 0, s, 0, c, 0, 0, 0, 0, 1];
}
function apply(m: Mat, p: number[]) {
  return [0, 1, 2].map((r) => m[r] * p[0] + m[4 + r] * p[1] + m[8 + r] * p[2] + m[12 + r]);
}
function lookAt(e: number[], t: number[], u: number[]): number[] {
  let zx = e[0] - t[0], zy = e[1] - t[1], zz = e[2] - t[2];
  let l = Math.hypot(zx, zy, zz); zx /= l; zy /= l; zz /= l;
  let xx = u[1] * zz - u[2] * zy, xy = u[2] * zx - u[0] * zz, xz = u[0] * zy - u[1] * zx;
  l = Math.hypot(xx, xy, xz); xx /= l; xy /= l; xz /= l;
  const yx = zy * xz - zz * xy, yy = zz * xx - zx * xz, yz = zx * xy - zy * xx;
  return [xx, yx, zx, 0, xy, yy, zy, 0, xz, yz, zz, 0,
    -(xx * e[0] + xy * e[1] + xz * e[2]), -(yx * e[0] + yy * e[1] + yz * e[2]), -(zx * e[0] + zy * e[1] + zz * e[2]), 1];
}
function mul(a: Mat, b: Mat) {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) {
    let s = 0;
    for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k];
    o[c * 4 + r] = s;
  }
  return o;
}

/** Unindexed positions + normals, smooth inside the crease angle only. */
export function shade(pos: Float32Array, idx: Uint32Array, creaseDeg = 40) {
  const T = idx.length / 3, V = pos.length / 3, p = pos;
  const fn = new Float64Array(3 * T), un = new Float64Array(3 * T);
  for (let t = 0; t < T; t++) {
    const a = 3 * idx[3 * t], b = 3 * idx[3 * t + 1], c = 3 * idx[3 * t + 2];
    const ux = p[b] - p[a], uy = p[b + 1] - p[a + 1], uz = p[b + 2] - p[a + 2];
    const vx = p[c] - p[a], vy = p[c + 1] - p[a + 1], vz = p[c + 2] - p[a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const l = Math.hypot(nx, ny, nz) || 1;
    fn[3 * t] = nx; fn[3 * t + 1] = ny; fn[3 * t + 2] = nz;
    un[3 * t] = nx / l; un[3 * t + 1] = ny / l; un[3 * t + 2] = nz / l;
  }
  const start = new Uint32Array(V + 1);
  for (let i = 0; i < idx.length; i++) start[idx[i] + 1]++;
  for (let v = 0; v < V; v++) start[v + 1] += start[v];
  const fill = start.slice(0, V), adj = new Uint32Array(idx.length);
  for (let i = 0; i < idx.length; i++) adj[fill[idx[i]]++] = (i / 3) | 0;
  const cosC = Math.cos((creaseDeg * Math.PI) / 180);
  const P = new Float32Array(9 * T), N = new Float32Array(9 * T);
  for (let t = 0; t < T; t++) for (let c = 0; c < 3; c++) {
    const v = idx[3 * t + c], o = 9 * t + 3 * c;
    let sx = 0, sy = 0, sz = 0;
    for (let k = start[v]; k < start[v + 1]; k++) {
      const s = adj[k];
      if (un[3 * t] * un[3 * s] + un[3 * t + 1] * un[3 * s + 1] + un[3 * t + 2] * un[3 * s + 2] > cosC) {
        sx += fn[3 * s]; sy += fn[3 * s + 1]; sz += fn[3 * s + 2];
      }
    }
    const l = Math.hypot(sx, sy, sz) || 1;
    P[o] = p[3 * v]; P[o + 1] = p[3 * v + 1]; P[o + 2] = p[3 * v + 2];
    N[o] = sx / l; N[o + 1] = sy / l; N[o + 2] = sz / l;
  }
  return { P, N };
}

function cssColor(name: string, fallback: number[]) {
  const s = getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  const m = /^#([0-9a-f]{6})$/i.exec(s);
  if (!m) return fallback;
  const v = parseInt(m[1], 16);
  return [((v >> 16) & 255) / 255, ((v >> 8) & 255) / 255, (v & 255) / 255];
}

export interface Scene {
  showWall: boolean;
  /** degrees past vertical: + overhang, - slab */
  wallAngle: number;
  ortho: boolean;
}

export interface Viewer {
  setMesh(pos: Float32Array, idx: Uint32Array): void;
  setScene(s: Partial<Scene>): void;
  redraw(): void;
  resetView(): void;
  /** look square at the wall face, or along it from the side */
  view(which: "front" | "side"): void;
}

const PMAX = 1.55;

function mix(a: number[], b: number[], t: number) {
  return a.map((v, i) => v + (b[i] - v) * t);
}

export function createViewer(canvas: HTMLCanvasElement): Viewer | null {
  let gl: WebGLRenderingContext | null = null;
  try { gl = canvas.getContext("webgl", { antialias: true, preserveDrawingBuffer: true }); } catch { gl = null; }
  if (!gl) return null;
  const g = gl;
  const prog = compile(g, VS, FS), lprog = compile(g, LVS, LFS);
  const loc = {
    aP: g.getAttribLocation(prog, "aP"), aN: g.getAttribLocation(prog, "aN"),
    uM: g.getUniformLocation(prog, "uM"), uCol: g.getUniformLocation(prog, "uCol"),
    uEye: g.getUniformLocation(prog, "uEye"), uKey: g.getUniformLocation(prog, "uKey"),
    uW: g.getUniformLocation(prog, "uW"),
    laP: g.getAttribLocation(lprog, "aP"), luM: g.getUniformLocation(lprog, "uM"), luCol: g.getUniformLocation(lprog, "uCol"),
  };
  const bufP = g.createBuffer(), bufN = g.createBuffer(), bufL = g.createBuffer();
  const bufWall = g.createBuffer(), bufGround = g.createBuffer();
  let count = 0, lineCount = 0;
  let centre = [0, 0, 0], radius = 60;
  let lo = [-1, -1, -1], hi = [1, 1, 1];
  let ext = { y: 60, z: 60 };
  const scene: Scene = { showWall: true, wallAngle: 0, ortho: false };
  // yaw 0 looks straight at the wall (camera on +x); the default is a
  // three-quarter view from the right and a little below the hold
  const home = { yaw: 0.65, pitch: 0.28 };
  const cam = { ...home, zoom: 1 };
  let pending = false;

  function setMesh(pos: Float32Array, idx: Uint32Array) {
    const { P, N } = shade(pos, idx);
    g.bindBuffer(g.ARRAY_BUFFER, bufP); g.bufferData(g.ARRAY_BUFFER, P, g.STATIC_DRAW);
    g.bindBuffer(g.ARRAY_BUFFER, bufN); g.bufferData(g.ARRAY_BUFFER, N, g.STATIC_DRAW);
    count = P.length / 3;
    lo = [Infinity, Infinity, Infinity]; hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < pos.length; i += 3) for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], pos[i + k]); hi[k] = Math.max(hi[k], pos[i + k]); }
    centre = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])];
    radius = 0.5 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
    // the wall: a 10 mm grid on x = 0 a little past the footprint
    const ey = Math.ceil((Math.max(Math.abs(lo[1]), Math.abs(hi[1])) + 20) / 10) * 10;
    const ez = Math.ceil((Math.max(Math.abs(lo[2]), Math.abs(hi[2])) + 20) / 10) * 10;
    ext = { y: ey, z: ez };
    // wall panel 1.5 mm behind the grid, so it hides the hold from behind (hold frame, turned by W in draw)
    g.bindBuffer(g.ARRAY_BUFFER, bufWall);
    g.bufferData(g.ARRAY_BUFFER, new Float32Array([
      -1.5, -ey, -ez, -1.5, ey, -ez, -1.5, ey, ez, -1.5, -ey, -ez, -1.5, ey, ez, -1.5, -ey, ez]), g.STATIC_DRAW);
    const gr: number[] = [];
    for (let y = -ey; y <= ey + 1e-6; y += 10) gr.push(-0.05, y, -ez, -0.05, y, ez);
    for (let z = -ez; z <= ez + 1e-6; z += 10) gr.push(-0.05, -ey, z, -0.05, ey, z);
    g.bindBuffer(g.ARRAY_BUFFER, bufL); g.bufferData(g.ARRAY_BUFFER, new Float32Array(gr), g.STATIC_DRAW);
    lineCount = gr.length / 3;
    redraw();
  }

  function draw() {
    pending = false;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = Math.max(1, Math.round(canvas.clientWidth * dpr)), h = Math.max(1, Math.round(canvas.clientHeight * dpr));
    if (canvas.width !== w || canvas.height !== h) { canvas.width = w; canvas.height = h; }
    g.viewport(0, 0, w, h);
    const bg = cssColor("--bg-soft", [0.95, 0.95, 0.96]);
    g.clearColor(bg[0], bg[1], bg[2], 1);
    g.clear(g.COLOR_BUFFER_BIT | g.DEPTH_BUFFER_BIT);
    g.enable(g.DEPTH_TEST);

    const fov = 0.6, asp = w / h;
    const W = rotY((scene.wallAngle * Math.PI) / 180);
    // the floor: horizontal, through the wall's bottom edge (or under the
    // hold, for a steep roof), reaching out past the hold
    let floor: number[] | null = null;
    let sceneR = radius;
    const t = apply(W, centre);
    if (scene.showWall) {
      const bot = apply(W, [-1.5, 0, -ext.z]);
      let zmin = Infinity;
      for (const x of [lo[0], hi[0]]) for (const z of [lo[2], hi[2]]) zmin = Math.min(zmin, apply(W, [x, 0, z])[2]);
      const fz = Math.min(bot[2], zmin - 10);
      const x0 = Math.min(bot[0], t[0] - radius), x1 = Math.max(bot[0], t[0]) + Math.max(2 * ext.z, 3 * radius);
      const y = ext.y;
      floor = [x0, -y, fz, x1, -y, fz, x1, y, fz, x0, -y, fz, x1, y, fz, x0, y, fz];
      for (let i = 0; i < 18; i += 3) sceneR = Math.max(sceneR, Math.hypot(floor[i] - t[0], floor[i + 1] - t[1], floor[i + 2] - t[2]));
      g.bindBuffer(g.ARRAY_BUFFER, bufGround); g.bufferData(g.ARRAY_BUFFER, new Float32Array(floor), g.DYNAMIC_DRAW);
    }
    sceneR = Math.max(sceneR, Math.hypot(ext.y, ext.z) + radius);

    const fit = radius / Math.sin(Math.min(fov, 2 * Math.atan(Math.tan(fov / 2) * asp)) / 2);
    const d = fit * 1.05 * cam.zoom;
    const cp = Math.cos(cam.pitch);
    const eye = [t[0] + d * cp * Math.cos(cam.yaw), t[1] + d * cp * Math.sin(cam.yaw), t[2] + d * Math.sin(cam.pitch)];
    // orthographic is sized to match the perspective view at the target
    const proj = scene.ortho
      ? ortho(d * Math.tan(fov / 2), asp, Math.max(0.5, d - 1.5 * sceneR), d + 1.5 * sceneR)
      : perspective(fov, asp, Math.max(0.5, d - 1.5 * sceneR), d + 1.5 * sceneR);
    const M = mul(proj, lookAt(eye, t, [0, 0, 1]));
    const MW = mul(M, W);

    const flat = (buf: WebGLBuffer | null, m: Float32Array, c: number[], mode: number, n: number) => {
      g.uniformMatrix4fv(loc.luM, false, m);
      g.uniform4f(loc.luCol, c[0], c[1], c[2], 1);
      g.bindBuffer(g.ARRAY_BUFFER, buf);
      g.enableVertexAttribArray(loc.laP);
      g.vertexAttribPointer(loc.laP, 3, g.FLOAT, false, 0, 0);
      g.drawArrays(mode, 0, n);
      g.disableVertexAttribArray(loc.laP);
    };
    const lc = cssColor("--line", [0.88, 0.89, 0.91]);
    const fg = cssColor("--fg", [0.1, 0.1, 0.12]);
    g.useProgram(lprog);
    if (scene.showWall && floor) {
      // push the panels back so the grid and the hold win the depth test
      g.enable(g.POLYGON_OFFSET_FILL); g.polygonOffset(1, 1);
      flat(bufWall, MW, mix(bg, fg, 0.07), g.TRIANGLES, 6);
      // from under the floor (front view of an overhang) it would hide everything
      if (eye[2] > floor[2]) flat(bufGround, M, mix(bg, fg, 0.15), g.TRIANGLES, 6);
      g.disable(g.POLYGON_OFFSET_FILL);
    }
    flat(bufL, MW, scene.showWall ? mix(bg, fg, 0.2) : lc, g.LINES, lineCount);

    if (!count) return;
    const col = cssColor("--accent", [0.18, 0.33, 0.78]).map((c) => 0.75 * c + 0.25);
    g.useProgram(prog);
    g.uniformMatrix4fv(loc.uM, false, M);
    g.uniformMatrix4fv(loc.uW, false, W);
    g.uniform3f(loc.uCol, col[0], col[1], col[2]);
    // in orthographic the view direction is fixed; put the "eye" far back
    const e2 = scene.ortho ? eye.map((v, i) => t[i] + (v - t[i]) * 50) : eye;
    g.uniform3f(loc.uEye, e2[0], e2[1], e2[2]);
    g.uniform3f(loc.uKey, 1.2, Math.sin(cam.yaw + 0.8), 1.4);
    g.bindBuffer(g.ARRAY_BUFFER, bufP); g.enableVertexAttribArray(loc.aP);
    g.vertexAttribPointer(loc.aP, 3, g.FLOAT, false, 0, 0);
    g.bindBuffer(g.ARRAY_BUFFER, bufN); g.enableVertexAttribArray(loc.aN);
    g.vertexAttribPointer(loc.aN, 3, g.FLOAT, false, 0, 0);
    g.drawArrays(g.TRIANGLES, 0, count);
    g.disableVertexAttribArray(loc.aP); g.disableVertexAttribArray(loc.aN);
  }
  function redraw() { if (!pending) { pending = true; requestAnimationFrame(draw); } }
  function resetView() { cam.yaw = home.yaw; cam.pitch = home.pitch; cam.zoom = 1; redraw(); }
  function view(which: "front" | "side") {
    if (which === "front") {
      // along the wall normal, W·(1,0,0) = (cos a, 0, -sin a)
      cam.yaw = 0;
      cam.pitch = Math.max(-PMAX, Math.min(PMAX, (-scene.wallAngle * Math.PI) / 180));
    } else {
      cam.yaw = -Math.PI / 2; cam.pitch = 0; // wall on the left, hold pointing right
    }
    redraw();
  }
  function setScene(s: Partial<Scene>) { Object.assign(scene, s); redraw(); }

  const pts = new Map<number, [number, number]>();
  let pinch0 = 0, zoom0 = 1;
  canvas.addEventListener("pointerdown", (e) => {
    canvas.setPointerCapture(e.pointerId);
    pts.set(e.pointerId, [e.clientX, e.clientY]);
    if (pts.size === 2) { const [a, b] = [...pts.values()]; pinch0 = Math.hypot(a[0] - b[0], a[1] - b[1]); zoom0 = cam.zoom; }
  });
  canvas.addEventListener("pointermove", (e) => {
    const prev = pts.get(e.pointerId);
    if (!prev) return;
    const cur: [number, number] = [e.clientX, e.clientY];
    pts.set(e.pointerId, cur);
    if (pts.size === 1) {
      cam.yaw -= (cur[0] - prev[0]) * 0.008;
      cam.pitch = Math.max(-PMAX, Math.min(PMAX, cam.pitch + (cur[1] - prev[1]) * 0.008));
    } else if (pts.size === 2 && pinch0 > 0) {
      const [a, b] = [...pts.values()];
      cam.zoom = Math.max(0.15, Math.min(4, (zoom0 * pinch0) / Math.max(1, Math.hypot(a[0] - b[0], a[1] - b[1]))));
    }
    redraw();
  });
  const up = (e: PointerEvent) => { pts.delete(e.pointerId); pinch0 = 0; };
  canvas.addEventListener("pointerup", up);
  canvas.addEventListener("pointercancel", up);
  canvas.addEventListener("wheel", (e) => {
    e.preventDefault();
    cam.zoom = Math.max(0.15, Math.min(4, cam.zoom * Math.exp(e.deltaY * 0.0015)));
    redraw();
  }, { passive: false });
  canvas.addEventListener("dblclick", resetView);
  new ResizeObserver(redraw).observe(canvas);
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", redraw);
  return { setMesh, setScene, redraw, resetView, view };
}
