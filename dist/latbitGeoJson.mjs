/**
 * @license
 * Copyright 2010-2023 Three.js Authors
 * SPDX-License-Identifier: MIT
 */
const ss = "155", Pn = { LEFT: 0, MIDDLE: 1, RIGHT: 2, ROTATE: 0, DOLLY: 1, PAN: 2 }, Ln = { ROTATE: 0, PAN: 1, DOLLY_PAN: 2, DOLLY_ROTATE: 3 }, Zo = 0, ys = 1, Ko = 2, $a = 1, Jo = 2, Yt = 3, ln = 0, gt = 1, Mt = 2, sn = 0, $n = 1, Ts = 2, As = 3, bs = 4, $o = 5, jn = 100, Qo = 101, el = 102, ws = 103, Rs = 104, tl = 200, nl = 201, il = 202, rl = 203, Qa = 204, eo = 205, sl = 206, al = 207, ol = 208, ll = 209, cl = 210, hl = 0, ul = 1, fl = 2, Yr = 3, dl = 4, pl = 5, ml = 6, gl = 7, to = 0, _l = 1, vl = 2, an = 0, xl = 1, Ml = 2, Sl = 3, El = 4, yl = 5, no = 300, ei = 301, ti = 302, qr = 303, jr = 304, rr = 306, Zr = 1e3, Lt = 1001, Kr = 1002, mt = 1003, Cs = 1004, dr = 1005, Tt = 1006, Tl = 1007, vi = 1008, on = 1009, Al = 1010, bl = 1011, as = 1012, io = 1013, nn = 1014, rn = 1015, xi = 1016, ro = 1017, so = 1018, xn = 1020, wl = 1021, Ut = 1023, Rl = 1024, Cl = 1025, Mn = 1026, ni = 1027, Pl = 1028, ao = 1029, Ll = 1030, oo = 1031, lo = 1033, pr = 33776, mr = 33777, gr = 33778, _r = 33779, Ps = 35840, Ls = 35841, Us = 35842, Ds = 35843, Ul = 36196, Is = 37492, Ns = 37496, Os = 37808, Fs = 37809, Bs = 37810, zs = 37811, Hs = 37812, Gs = 37813, Vs = 37814, ks = 37815, Ws = 37816, Xs = 37817, Ys = 37818, qs = 37819, js = 37820, Zs = 37821, vr = 36492, Dl = 36283, Ks = 36284, Js = 36285, $s = 36286, co = 3e3, Sn = 3001, Il = 3200, Nl = 3201, ho = 0, Ol = 1, En = "", Oe = "srgb", Ft = "srgb-linear", uo = "display-p3", xr = 7680, Fl = 519, Bl = 512, zl = 513, Hl = 514, Gl = 515, Vl = 516, kl = 517, Wl = 518, Xl = 519, Qs = 35044, ea = "300 es", Jr = 1035, qt = 2e3, tr = 2001;
class Rn {
  addEventListener(e, t) {
    this._listeners === void 0 && (this._listeners = {});
    const n = this._listeners;
    n[e] === void 0 && (n[e] = []), n[e].indexOf(t) === -1 && n[e].push(t);
  }
  hasEventListener(e, t) {
    if (this._listeners === void 0)
      return !1;
    const n = this._listeners;
    return n[e] !== void 0 && n[e].indexOf(t) !== -1;
  }
  removeEventListener(e, t) {
    if (this._listeners === void 0)
      return;
    const r = this._listeners[e];
    if (r !== void 0) {
      const s = r.indexOf(t);
      s !== -1 && r.splice(s, 1);
    }
  }
  dispatchEvent(e) {
    if (this._listeners === void 0)
      return;
    const n = this._listeners[e.type];
    if (n !== void 0) {
      e.target = this;
      const r = n.slice(0);
      for (let s = 0, o = r.length; s < o; s++)
        r[s].call(this, e);
      e.target = null;
    }
  }
}
const lt = ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "0a", "0b", "0c", "0d", "0e", "0f", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "1a", "1b", "1c", "1d", "1e", "1f", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "2a", "2b", "2c", "2d", "2e", "2f", "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "3a", "3b", "3c", "3d", "3e", "3f", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "4a", "4b", "4c", "4d", "4e", "4f", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59", "5a", "5b", "5c", "5d", "5e", "5f", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "6a", "6b", "6c", "6d", "6e", "6f", "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "7a", "7b", "7c", "7d", "7e", "7f", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "8a", "8b", "8c", "8d", "8e", "8f", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "9a", "9b", "9c", "9d", "9e", "9f", "a0", "a1", "a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "aa", "ab", "ac", "ad", "ae", "af", "b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "ba", "bb", "bc", "bd", "be", "bf", "c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9", "ca", "cb", "cc", "cd", "ce", "cf", "d0", "d1", "d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "da", "db", "dc", "dd", "de", "df", "e0", "e1", "e2", "e3", "e4", "e5", "e6", "e7", "e8", "e9", "ea", "eb", "ec", "ed", "ee", "ef", "f0", "f1", "f2", "f3", "f4", "f5", "f6", "f7", "f8", "f9", "fa", "fb", "fc", "fd", "fe", "ff"];
let ta = 1234567;
const ui = Math.PI / 180, Mi = 180 / Math.PI;
function Cn() {
  const i = Math.random() * 4294967295 | 0, e = Math.random() * 4294967295 | 0, t = Math.random() * 4294967295 | 0, n = Math.random() * 4294967295 | 0;
  return (lt[i & 255] + lt[i >> 8 & 255] + lt[i >> 16 & 255] + lt[i >> 24 & 255] + "-" + lt[e & 255] + lt[e >> 8 & 255] + "-" + lt[e >> 16 & 15 | 64] + lt[e >> 24 & 255] + "-" + lt[t & 63 | 128] + lt[t >> 8 & 255] + "-" + lt[t >> 16 & 255] + lt[t >> 24 & 255] + lt[n & 255] + lt[n >> 8 & 255] + lt[n >> 16 & 255] + lt[n >> 24 & 255]).toLowerCase();
}
function it(i, e, t) {
  return Math.max(e, Math.min(t, i));
}
function os(i, e) {
  return (i % e + e) % e;
}
function Yl(i, e, t, n, r) {
  return n + (i - e) * (r - n) / (t - e);
}
function ql(i, e, t) {
  return i !== e ? (t - i) / (e - i) : 0;
}
function fi(i, e, t) {
  return (1 - t) * i + t * e;
}
function jl(i, e, t, n) {
  return fi(i, e, 1 - Math.exp(-t * n));
}
function Zl(i, e = 1) {
  return e - Math.abs(os(i, e * 2) - e);
}
function Kl(i, e, t) {
  return i <= e ? 0 : i >= t ? 1 : (i = (i - e) / (t - e), i * i * (3 - 2 * i));
}
function Jl(i, e, t) {
  return i <= e ? 0 : i >= t ? 1 : (i = (i - e) / (t - e), i * i * i * (i * (i * 6 - 15) + 10));
}
function $l(i, e) {
  return i + Math.floor(Math.random() * (e - i + 1));
}
function Ql(i, e) {
  return i + Math.random() * (e - i);
}
function ec(i) {
  return i * (0.5 - Math.random());
}
function tc(i) {
  i !== void 0 && (ta = i);
  let e = ta += 1831565813;
  return e = Math.imul(e ^ e >>> 15, e | 1), e ^= e + Math.imul(e ^ e >>> 7, e | 61), ((e ^ e >>> 14) >>> 0) / 4294967296;
}
function nc(i) {
  return i * ui;
}
function ic(i) {
  return i * Mi;
}
function $r(i) {
  return (i & i - 1) === 0 && i !== 0;
}
function rc(i) {
  return Math.pow(2, Math.ceil(Math.log(i) / Math.LN2));
}
function nr(i) {
  return Math.pow(2, Math.floor(Math.log(i) / Math.LN2));
}
function sc(i, e, t, n, r) {
  const s = Math.cos, o = Math.sin, a = s(t / 2), l = o(t / 2), c = s((e + n) / 2), h = o((e + n) / 2), f = s((e - n) / 2), u = o((e - n) / 2), m = s((n - e) / 2), g = o((n - e) / 2);
  switch (r) {
    case "XYX":
      i.set(a * h, l * f, l * u, a * c);
      break;
    case "YZY":
      i.set(l * u, a * h, l * f, a * c);
      break;
    case "ZXZ":
      i.set(l * f, l * u, a * h, a * c);
      break;
    case "XZX":
      i.set(a * h, l * g, l * m, a * c);
      break;
    case "YXY":
      i.set(l * m, a * h, l * g, a * c);
      break;
    case "ZYZ":
      i.set(l * g, l * m, a * h, a * c);
      break;
    default:
      console.warn("THREE.MathUtils: .setQuaternionFromProperEuler() encountered an unknown order: " + r);
  }
}
function Zn(i, e) {
  switch (e.constructor) {
    case Float32Array:
      return i;
    case Uint32Array:
      return i / 4294967295;
    case Uint16Array:
      return i / 65535;
    case Uint8Array:
      return i / 255;
    case Int32Array:
      return Math.max(i / 2147483647, -1);
    case Int16Array:
      return Math.max(i / 32767, -1);
    case Int8Array:
      return Math.max(i / 127, -1);
    default:
      throw new Error("Invalid component type.");
  }
}
function dt(i, e) {
  switch (e.constructor) {
    case Float32Array:
      return i;
    case Uint32Array:
      return Math.round(i * 4294967295);
    case Uint16Array:
      return Math.round(i * 65535);
    case Uint8Array:
      return Math.round(i * 255);
    case Int32Array:
      return Math.round(i * 2147483647);
    case Int16Array:
      return Math.round(i * 32767);
    case Int8Array:
      return Math.round(i * 127);
    default:
      throw new Error("Invalid component type.");
  }
}
const ac = {
  DEG2RAD: ui,
  RAD2DEG: Mi,
  generateUUID: Cn,
  clamp: it,
  euclideanModulo: os,
  mapLinear: Yl,
  inverseLerp: ql,
  lerp: fi,
  damp: jl,
  pingpong: Zl,
  smoothstep: Kl,
  smootherstep: Jl,
  randInt: $l,
  randFloat: Ql,
  randFloatSpread: ec,
  seededRandom: tc,
  degToRad: nc,
  radToDeg: ic,
  isPowerOfTwo: $r,
  ceilPowerOfTwo: rc,
  floorPowerOfTwo: nr,
  setQuaternionFromProperEuler: sc,
  normalize: dt,
  denormalize: Zn
};
class oe {
  constructor(e = 0, t = 0) {
    oe.prototype.isVector2 = !0, this.x = e, this.y = t;
  }
  get width() {
    return this.x;
  }
  set width(e) {
    this.x = e;
  }
  get height() {
    return this.y;
  }
  set height(e) {
    this.y = e;
  }
  set(e, t) {
    return this.x = e, this.y = t, this;
  }
  setScalar(e) {
    return this.x = e, this.y = e, this;
  }
  setX(e) {
    return this.x = e, this;
  }
  setY(e) {
    return this.y = e, this;
  }
  setComponent(e, t) {
    switch (e) {
      case 0:
        this.x = t;
        break;
      case 1:
        this.y = t;
        break;
      default:
        throw new Error("index is out of range: " + e);
    }
    return this;
  }
  getComponent(e) {
    switch (e) {
      case 0:
        return this.x;
      case 1:
        return this.y;
      default:
        throw new Error("index is out of range: " + e);
    }
  }
  clone() {
    return new this.constructor(this.x, this.y);
  }
  copy(e) {
    return this.x = e.x, this.y = e.y, this;
  }
  add(e) {
    return this.x += e.x, this.y += e.y, this;
  }
  addScalar(e) {
    return this.x += e, this.y += e, this;
  }
  addVectors(e, t) {
    return this.x = e.x + t.x, this.y = e.y + t.y, this;
  }
  addScaledVector(e, t) {
    return this.x += e.x * t, this.y += e.y * t, this;
  }
  sub(e) {
    return this.x -= e.x, this.y -= e.y, this;
  }
  subScalar(e) {
    return this.x -= e, this.y -= e, this;
  }
  subVectors(e, t) {
    return this.x = e.x - t.x, this.y = e.y - t.y, this;
  }
  multiply(e) {
    return this.x *= e.x, this.y *= e.y, this;
  }
  multiplyScalar(e) {
    return this.x *= e, this.y *= e, this;
  }
  divide(e) {
    return this.x /= e.x, this.y /= e.y, this;
  }
  divideScalar(e) {
    return this.multiplyScalar(1 / e);
  }
  applyMatrix3(e) {
    const t = this.x, n = this.y, r = e.elements;
    return this.x = r[0] * t + r[3] * n + r[6], this.y = r[1] * t + r[4] * n + r[7], this;
  }
  min(e) {
    return this.x = Math.min(this.x, e.x), this.y = Math.min(this.y, e.y), this;
  }
  max(e) {
    return this.x = Math.max(this.x, e.x), this.y = Math.max(this.y, e.y), this;
  }
  clamp(e, t) {
    return this.x = Math.max(e.x, Math.min(t.x, this.x)), this.y = Math.max(e.y, Math.min(t.y, this.y)), this;
  }
  clampScalar(e, t) {
    return this.x = Math.max(e, Math.min(t, this.x)), this.y = Math.max(e, Math.min(t, this.y)), this;
  }
  clampLength(e, t) {
    const n = this.length();
    return this.divideScalar(n || 1).multiplyScalar(Math.max(e, Math.min(t, n)));
  }
  floor() {
    return this.x = Math.floor(this.x), this.y = Math.floor(this.y), this;
  }
  ceil() {
    return this.x = Math.ceil(this.x), this.y = Math.ceil(this.y), this;
  }
  round() {
    return this.x = Math.round(this.x), this.y = Math.round(this.y), this;
  }
  roundToZero() {
    return this.x = this.x < 0 ? Math.ceil(this.x) : Math.floor(this.x), this.y = this.y < 0 ? Math.ceil(this.y) : Math.floor(this.y), this;
  }
  negate() {
    return this.x = -this.x, this.y = -this.y, this;
  }
  dot(e) {
    return this.x * e.x + this.y * e.y;
  }
  cross(e) {
    return this.x * e.y - this.y * e.x;
  }
  lengthSq() {
    return this.x * this.x + this.y * this.y;
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y);
  }
  manhattanLength() {
    return Math.abs(this.x) + Math.abs(this.y);
  }
  normalize() {
    return this.divideScalar(this.length() || 1);
  }
  angle() {
    return Math.atan2(-this.y, -this.x) + Math.PI;
  }
  angleTo(e) {
    const t = Math.sqrt(this.lengthSq() * e.lengthSq());
    if (t === 0)
      return Math.PI / 2;
    const n = this.dot(e) / t;
    return Math.acos(it(n, -1, 1));
  }
  distanceTo(e) {
    return Math.sqrt(this.distanceToSquared(e));
  }
  distanceToSquared(e) {
    const t = this.x - e.x, n = this.y - e.y;
    return t * t + n * n;
  }
  manhattanDistanceTo(e) {
    return Math.abs(this.x - e.x) + Math.abs(this.y - e.y);
  }
  setLength(e) {
    return this.normalize().multiplyScalar(e);
  }
  lerp(e, t) {
    return this.x += (e.x - this.x) * t, this.y += (e.y - this.y) * t, this;
  }
  lerpVectors(e, t, n) {
    return this.x = e.x + (t.x - e.x) * n, this.y = e.y + (t.y - e.y) * n, this;
  }
  equals(e) {
    return e.x === this.x && e.y === this.y;
  }
  fromArray(e, t = 0) {
    return this.x = e[t], this.y = e[t + 1], this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this.x, e[t + 1] = this.y, e;
  }
  fromBufferAttribute(e, t) {
    return this.x = e.getX(t), this.y = e.getY(t), this;
  }
  rotateAround(e, t) {
    const n = Math.cos(t), r = Math.sin(t), s = this.x - e.x, o = this.y - e.y;
    return this.x = s * n - o * r + e.x, this.y = s * r + o * n + e.y, this;
  }
  random() {
    return this.x = Math.random(), this.y = Math.random(), this;
  }
  *[Symbol.iterator]() {
    yield this.x, yield this.y;
  }
}
class Be {
  constructor(e, t, n, r, s, o, a, l, c) {
    Be.prototype.isMatrix3 = !0, this.elements = [
      1,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      1
    ], e !== void 0 && this.set(e, t, n, r, s, o, a, l, c);
  }
  set(e, t, n, r, s, o, a, l, c) {
    const h = this.elements;
    return h[0] = e, h[1] = r, h[2] = a, h[3] = t, h[4] = s, h[5] = l, h[6] = n, h[7] = o, h[8] = c, this;
  }
  identity() {
    return this.set(
      1,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      1
    ), this;
  }
  copy(e) {
    const t = this.elements, n = e.elements;
    return t[0] = n[0], t[1] = n[1], t[2] = n[2], t[3] = n[3], t[4] = n[4], t[5] = n[5], t[6] = n[6], t[7] = n[7], t[8] = n[8], this;
  }
  extractBasis(e, t, n) {
    return e.setFromMatrix3Column(this, 0), t.setFromMatrix3Column(this, 1), n.setFromMatrix3Column(this, 2), this;
  }
  setFromMatrix4(e) {
    const t = e.elements;
    return this.set(
      t[0],
      t[4],
      t[8],
      t[1],
      t[5],
      t[9],
      t[2],
      t[6],
      t[10]
    ), this;
  }
  multiply(e) {
    return this.multiplyMatrices(this, e);
  }
  premultiply(e) {
    return this.multiplyMatrices(e, this);
  }
  multiplyMatrices(e, t) {
    const n = e.elements, r = t.elements, s = this.elements, o = n[0], a = n[3], l = n[6], c = n[1], h = n[4], f = n[7], u = n[2], m = n[5], g = n[8], x = r[0], p = r[3], d = r[6], A = r[1], _ = r[4], T = r[7], C = r[2], L = r[5], w = r[8];
    return s[0] = o * x + a * A + l * C, s[3] = o * p + a * _ + l * L, s[6] = o * d + a * T + l * w, s[1] = c * x + h * A + f * C, s[4] = c * p + h * _ + f * L, s[7] = c * d + h * T + f * w, s[2] = u * x + m * A + g * C, s[5] = u * p + m * _ + g * L, s[8] = u * d + m * T + g * w, this;
  }
  multiplyScalar(e) {
    const t = this.elements;
    return t[0] *= e, t[3] *= e, t[6] *= e, t[1] *= e, t[4] *= e, t[7] *= e, t[2] *= e, t[5] *= e, t[8] *= e, this;
  }
  determinant() {
    const e = this.elements, t = e[0], n = e[1], r = e[2], s = e[3], o = e[4], a = e[5], l = e[6], c = e[7], h = e[8];
    return t * o * h - t * a * c - n * s * h + n * a * l + r * s * c - r * o * l;
  }
  invert() {
    const e = this.elements, t = e[0], n = e[1], r = e[2], s = e[3], o = e[4], a = e[5], l = e[6], c = e[7], h = e[8], f = h * o - a * c, u = a * l - h * s, m = c * s - o * l, g = t * f + n * u + r * m;
    if (g === 0)
      return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0);
    const x = 1 / g;
    return e[0] = f * x, e[1] = (r * c - h * n) * x, e[2] = (a * n - r * o) * x, e[3] = u * x, e[4] = (h * t - r * l) * x, e[5] = (r * s - a * t) * x, e[6] = m * x, e[7] = (n * l - c * t) * x, e[8] = (o * t - n * s) * x, this;
  }
  transpose() {
    let e;
    const t = this.elements;
    return e = t[1], t[1] = t[3], t[3] = e, e = t[2], t[2] = t[6], t[6] = e, e = t[5], t[5] = t[7], t[7] = e, this;
  }
  getNormalMatrix(e) {
    return this.setFromMatrix4(e).invert().transpose();
  }
  transposeIntoArray(e) {
    const t = this.elements;
    return e[0] = t[0], e[1] = t[3], e[2] = t[6], e[3] = t[1], e[4] = t[4], e[5] = t[7], e[6] = t[2], e[7] = t[5], e[8] = t[8], this;
  }
  setUvTransform(e, t, n, r, s, o, a) {
    const l = Math.cos(s), c = Math.sin(s);
    return this.set(
      n * l,
      n * c,
      -n * (l * o + c * a) + o + e,
      -r * c,
      r * l,
      -r * (-c * o + l * a) + a + t,
      0,
      0,
      1
    ), this;
  }
  //
  scale(e, t) {
    return this.premultiply(Mr.makeScale(e, t)), this;
  }
  rotate(e) {
    return this.premultiply(Mr.makeRotation(-e)), this;
  }
  translate(e, t) {
    return this.premultiply(Mr.makeTranslation(e, t)), this;
  }
  // for 2D Transforms
  makeTranslation(e, t) {
    return e.isVector2 ? this.set(
      1,
      0,
      e.x,
      0,
      1,
      e.y,
      0,
      0,
      1
    ) : this.set(
      1,
      0,
      e,
      0,
      1,
      t,
      0,
      0,
      1
    ), this;
  }
  makeRotation(e) {
    const t = Math.cos(e), n = Math.sin(e);
    return this.set(
      t,
      -n,
      0,
      n,
      t,
      0,
      0,
      0,
      1
    ), this;
  }
  makeScale(e, t) {
    return this.set(
      e,
      0,
      0,
      0,
      t,
      0,
      0,
      0,
      1
    ), this;
  }
  //
  equals(e) {
    const t = this.elements, n = e.elements;
    for (let r = 0; r < 9; r++)
      if (t[r] !== n[r])
        return !1;
    return !0;
  }
  fromArray(e, t = 0) {
    for (let n = 0; n < 9; n++)
      this.elements[n] = e[n + t];
    return this;
  }
  toArray(e = [], t = 0) {
    const n = this.elements;
    return e[t] = n[0], e[t + 1] = n[1], e[t + 2] = n[2], e[t + 3] = n[3], e[t + 4] = n[4], e[t + 5] = n[5], e[t + 6] = n[6], e[t + 7] = n[7], e[t + 8] = n[8], e;
  }
  clone() {
    return new this.constructor().fromArray(this.elements);
  }
}
const Mr = /* @__PURE__ */ new Be();
function fo(i) {
  for (let e = i.length - 1; e >= 0; --e)
    if (i[e] >= 65535)
      return !0;
  return !1;
}
function ir(i) {
  return document.createElementNS("http://www.w3.org/1999/xhtml", i);
}
const na = {};
function di(i) {
  i in na || (na[i] = !0, console.warn(i));
}
function Qn(i) {
  return i < 0.04045 ? i * 0.0773993808 : Math.pow(i * 0.9478672986 + 0.0521327014, 2.4);
}
function Sr(i) {
  return i < 31308e-7 ? i * 12.92 : 1.055 * Math.pow(i, 0.41666) - 0.055;
}
const oc = /* @__PURE__ */ new Be().fromArray([
  0.8224621,
  0.0331941,
  0.0170827,
  0.177538,
  0.9668058,
  0.0723974,
  -1e-7,
  1e-7,
  0.9105199
]), lc = /* @__PURE__ */ new Be().fromArray([
  1.2249401,
  -0.0420569,
  -0.0196376,
  -0.2249404,
  1.0420571,
  -0.0786361,
  1e-7,
  0,
  1.0982735
]);
function cc(i) {
  return i.convertSRGBToLinear().applyMatrix3(lc);
}
function hc(i) {
  return i.applyMatrix3(oc).convertLinearToSRGB();
}
const uc = {
  [Ft]: (i) => i,
  [Oe]: (i) => i.convertSRGBToLinear(),
  [uo]: cc
}, fc = {
  [Ft]: (i) => i,
  [Oe]: (i) => i.convertLinearToSRGB(),
  [uo]: hc
}, bt = {
  enabled: !0,
  get legacyMode() {
    return console.warn("THREE.ColorManagement: .legacyMode=false renamed to .enabled=true in r150."), !this.enabled;
  },
  set legacyMode(i) {
    console.warn("THREE.ColorManagement: .legacyMode=false renamed to .enabled=true in r150."), this.enabled = !i;
  },
  get workingColorSpace() {
    return Ft;
  },
  set workingColorSpace(i) {
    console.warn("THREE.ColorManagement: .workingColorSpace is readonly.");
  },
  convert: function(i, e, t) {
    if (this.enabled === !1 || e === t || !e || !t)
      return i;
    const n = uc[e], r = fc[t];
    if (n === void 0 || r === void 0)
      throw new Error(`Unsupported color space conversion, "${e}" to "${t}".`);
    return r(n(i));
  },
  fromWorkingColorSpace: function(i, e) {
    return this.convert(i, this.workingColorSpace, e);
  },
  toWorkingColorSpace: function(i, e) {
    return this.convert(i, e, this.workingColorSpace);
  }
};
let Un;
class po {
  static getDataURL(e) {
    if (/^data:/i.test(e.src) || typeof HTMLCanvasElement > "u")
      return e.src;
    let t;
    if (e instanceof HTMLCanvasElement)
      t = e;
    else {
      Un === void 0 && (Un = ir("canvas")), Un.width = e.width, Un.height = e.height;
      const n = Un.getContext("2d");
      e instanceof ImageData ? n.putImageData(e, 0, 0) : n.drawImage(e, 0, 0, e.width, e.height), t = Un;
    }
    return t.width > 2048 || t.height > 2048 ? (console.warn("THREE.ImageUtils.getDataURL: Image converted to jpg for performance reasons", e), t.toDataURL("image/jpeg", 0.6)) : t.toDataURL("image/png");
  }
  static sRGBToLinear(e) {
    if (typeof HTMLImageElement < "u" && e instanceof HTMLImageElement || typeof HTMLCanvasElement < "u" && e instanceof HTMLCanvasElement || typeof ImageBitmap < "u" && e instanceof ImageBitmap) {
      const t = ir("canvas");
      t.width = e.width, t.height = e.height;
      const n = t.getContext("2d");
      n.drawImage(e, 0, 0, e.width, e.height);
      const r = n.getImageData(0, 0, e.width, e.height), s = r.data;
      for (let o = 0; o < s.length; o++)
        s[o] = Qn(s[o] / 255) * 255;
      return n.putImageData(r, 0, 0), t;
    } else if (e.data) {
      const t = e.data.slice(0);
      for (let n = 0; n < t.length; n++)
        t instanceof Uint8Array || t instanceof Uint8ClampedArray ? t[n] = Math.floor(Qn(t[n] / 255) * 255) : t[n] = Qn(t[n]);
      return {
        data: t,
        width: e.width,
        height: e.height
      };
    } else
      return console.warn("THREE.ImageUtils.sRGBToLinear(): Unsupported image type. No color space conversion applied."), e;
  }
}
let dc = 0;
class mo {
  constructor(e = null) {
    this.isSource = !0, Object.defineProperty(this, "id", { value: dc++ }), this.uuid = Cn(), this.data = e, this.version = 0;
  }
  set needsUpdate(e) {
    e === !0 && this.version++;
  }
  toJSON(e) {
    const t = e === void 0 || typeof e == "string";
    if (!t && e.images[this.uuid] !== void 0)
      return e.images[this.uuid];
    const n = {
      uuid: this.uuid,
      url: ""
    }, r = this.data;
    if (r !== null) {
      let s;
      if (Array.isArray(r)) {
        s = [];
        for (let o = 0, a = r.length; o < a; o++)
          r[o].isDataTexture ? s.push(Er(r[o].image)) : s.push(Er(r[o]));
      } else
        s = Er(r);
      n.url = s;
    }
    return t || (e.images[this.uuid] = n), n;
  }
}
function Er(i) {
  return typeof HTMLImageElement < "u" && i instanceof HTMLImageElement || typeof HTMLCanvasElement < "u" && i instanceof HTMLCanvasElement || typeof ImageBitmap < "u" && i instanceof ImageBitmap ? po.getDataURL(i) : i.data ? {
    data: Array.from(i.data),
    width: i.width,
    height: i.height,
    type: i.data.constructor.name
  } : (console.warn("THREE.Texture: Unable to serialize Texture."), {});
}
let pc = 0;
class St extends Rn {
  constructor(e = St.DEFAULT_IMAGE, t = St.DEFAULT_MAPPING, n = Lt, r = Lt, s = Tt, o = vi, a = Ut, l = on, c = St.DEFAULT_ANISOTROPY, h = En) {
    super(), this.isTexture = !0, Object.defineProperty(this, "id", { value: pc++ }), this.uuid = Cn(), this.name = "", this.source = new mo(e), this.mipmaps = [], this.mapping = t, this.channel = 0, this.wrapS = n, this.wrapT = r, this.magFilter = s, this.minFilter = o, this.anisotropy = c, this.format = a, this.internalFormat = null, this.type = l, this.offset = new oe(0, 0), this.repeat = new oe(1, 1), this.center = new oe(0, 0), this.rotation = 0, this.matrixAutoUpdate = !0, this.matrix = new Be(), this.generateMipmaps = !0, this.premultiplyAlpha = !1, this.flipY = !0, this.unpackAlignment = 4, typeof h == "string" ? this.colorSpace = h : (di("THREE.Texture: Property .encoding has been replaced by .colorSpace."), this.colorSpace = h === Sn ? Oe : En), this.userData = {}, this.version = 0, this.onUpdate = null, this.isRenderTargetTexture = !1, this.needsPMREMUpdate = !1;
  }
  get image() {
    return this.source.data;
  }
  set image(e = null) {
    this.source.data = e;
  }
  updateMatrix() {
    this.matrix.setUvTransform(this.offset.x, this.offset.y, this.repeat.x, this.repeat.y, this.rotation, this.center.x, this.center.y);
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    return this.name = e.name, this.source = e.source, this.mipmaps = e.mipmaps.slice(0), this.mapping = e.mapping, this.channel = e.channel, this.wrapS = e.wrapS, this.wrapT = e.wrapT, this.magFilter = e.magFilter, this.minFilter = e.minFilter, this.anisotropy = e.anisotropy, this.format = e.format, this.internalFormat = e.internalFormat, this.type = e.type, this.offset.copy(e.offset), this.repeat.copy(e.repeat), this.center.copy(e.center), this.rotation = e.rotation, this.matrixAutoUpdate = e.matrixAutoUpdate, this.matrix.copy(e.matrix), this.generateMipmaps = e.generateMipmaps, this.premultiplyAlpha = e.premultiplyAlpha, this.flipY = e.flipY, this.unpackAlignment = e.unpackAlignment, this.colorSpace = e.colorSpace, this.userData = JSON.parse(JSON.stringify(e.userData)), this.needsUpdate = !0, this;
  }
  toJSON(e) {
    const t = e === void 0 || typeof e == "string";
    if (!t && e.textures[this.uuid] !== void 0)
      return e.textures[this.uuid];
    const n = {
      metadata: {
        version: 4.6,
        type: "Texture",
        generator: "Texture.toJSON"
      },
      uuid: this.uuid,
      name: this.name,
      image: this.source.toJSON(e).uuid,
      mapping: this.mapping,
      channel: this.channel,
      repeat: [this.repeat.x, this.repeat.y],
      offset: [this.offset.x, this.offset.y],
      center: [this.center.x, this.center.y],
      rotation: this.rotation,
      wrap: [this.wrapS, this.wrapT],
      format: this.format,
      internalFormat: this.internalFormat,
      type: this.type,
      colorSpace: this.colorSpace,
      minFilter: this.minFilter,
      magFilter: this.magFilter,
      anisotropy: this.anisotropy,
      flipY: this.flipY,
      generateMipmaps: this.generateMipmaps,
      premultiplyAlpha: this.premultiplyAlpha,
      unpackAlignment: this.unpackAlignment
    };
    return Object.keys(this.userData).length > 0 && (n.userData = this.userData), t || (e.textures[this.uuid] = n), n;
  }
  dispose() {
    this.dispatchEvent({ type: "dispose" });
  }
  transformUv(e) {
    if (this.mapping !== no)
      return e;
    if (e.applyMatrix3(this.matrix), e.x < 0 || e.x > 1)
      switch (this.wrapS) {
        case Zr:
          e.x = e.x - Math.floor(e.x);
          break;
        case Lt:
          e.x = e.x < 0 ? 0 : 1;
          break;
        case Kr:
          Math.abs(Math.floor(e.x) % 2) === 1 ? e.x = Math.ceil(e.x) - e.x : e.x = e.x - Math.floor(e.x);
          break;
      }
    if (e.y < 0 || e.y > 1)
      switch (this.wrapT) {
        case Zr:
          e.y = e.y - Math.floor(e.y);
          break;
        case Lt:
          e.y = e.y < 0 ? 0 : 1;
          break;
        case Kr:
          Math.abs(Math.floor(e.y) % 2) === 1 ? e.y = Math.ceil(e.y) - e.y : e.y = e.y - Math.floor(e.y);
          break;
      }
    return this.flipY && (e.y = 1 - e.y), e;
  }
  set needsUpdate(e) {
    e === !0 && (this.version++, this.source.needsUpdate = !0);
  }
  get encoding() {
    return di("THREE.Texture: Property .encoding has been replaced by .colorSpace."), this.colorSpace === Oe ? Sn : co;
  }
  set encoding(e) {
    di("THREE.Texture: Property .encoding has been replaced by .colorSpace."), this.colorSpace = e === Sn ? Oe : En;
  }
}
St.DEFAULT_IMAGE = null;
St.DEFAULT_MAPPING = no;
St.DEFAULT_ANISOTROPY = 1;
class ot {
  constructor(e = 0, t = 0, n = 0, r = 1) {
    ot.prototype.isVector4 = !0, this.x = e, this.y = t, this.z = n, this.w = r;
  }
  get width() {
    return this.z;
  }
  set width(e) {
    this.z = e;
  }
  get height() {
    return this.w;
  }
  set height(e) {
    this.w = e;
  }
  set(e, t, n, r) {
    return this.x = e, this.y = t, this.z = n, this.w = r, this;
  }
  setScalar(e) {
    return this.x = e, this.y = e, this.z = e, this.w = e, this;
  }
  setX(e) {
    return this.x = e, this;
  }
  setY(e) {
    return this.y = e, this;
  }
  setZ(e) {
    return this.z = e, this;
  }
  setW(e) {
    return this.w = e, this;
  }
  setComponent(e, t) {
    switch (e) {
      case 0:
        this.x = t;
        break;
      case 1:
        this.y = t;
        break;
      case 2:
        this.z = t;
        break;
      case 3:
        this.w = t;
        break;
      default:
        throw new Error("index is out of range: " + e);
    }
    return this;
  }
  getComponent(e) {
    switch (e) {
      case 0:
        return this.x;
      case 1:
        return this.y;
      case 2:
        return this.z;
      case 3:
        return this.w;
      default:
        throw new Error("index is out of range: " + e);
    }
  }
  clone() {
    return new this.constructor(this.x, this.y, this.z, this.w);
  }
  copy(e) {
    return this.x = e.x, this.y = e.y, this.z = e.z, this.w = e.w !== void 0 ? e.w : 1, this;
  }
  add(e) {
    return this.x += e.x, this.y += e.y, this.z += e.z, this.w += e.w, this;
  }
  addScalar(e) {
    return this.x += e, this.y += e, this.z += e, this.w += e, this;
  }
  addVectors(e, t) {
    return this.x = e.x + t.x, this.y = e.y + t.y, this.z = e.z + t.z, this.w = e.w + t.w, this;
  }
  addScaledVector(e, t) {
    return this.x += e.x * t, this.y += e.y * t, this.z += e.z * t, this.w += e.w * t, this;
  }
  sub(e) {
    return this.x -= e.x, this.y -= e.y, this.z -= e.z, this.w -= e.w, this;
  }
  subScalar(e) {
    return this.x -= e, this.y -= e, this.z -= e, this.w -= e, this;
  }
  subVectors(e, t) {
    return this.x = e.x - t.x, this.y = e.y - t.y, this.z = e.z - t.z, this.w = e.w - t.w, this;
  }
  multiply(e) {
    return this.x *= e.x, this.y *= e.y, this.z *= e.z, this.w *= e.w, this;
  }
  multiplyScalar(e) {
    return this.x *= e, this.y *= e, this.z *= e, this.w *= e, this;
  }
  applyMatrix4(e) {
    const t = this.x, n = this.y, r = this.z, s = this.w, o = e.elements;
    return this.x = o[0] * t + o[4] * n + o[8] * r + o[12] * s, this.y = o[1] * t + o[5] * n + o[9] * r + o[13] * s, this.z = o[2] * t + o[6] * n + o[10] * r + o[14] * s, this.w = o[3] * t + o[7] * n + o[11] * r + o[15] * s, this;
  }
  divideScalar(e) {
    return this.multiplyScalar(1 / e);
  }
  setAxisAngleFromQuaternion(e) {
    this.w = 2 * Math.acos(e.w);
    const t = Math.sqrt(1 - e.w * e.w);
    return t < 1e-4 ? (this.x = 1, this.y = 0, this.z = 0) : (this.x = e.x / t, this.y = e.y / t, this.z = e.z / t), this;
  }
  setAxisAngleFromRotationMatrix(e) {
    let t, n, r, s;
    const l = e.elements, c = l[0], h = l[4], f = l[8], u = l[1], m = l[5], g = l[9], x = l[2], p = l[6], d = l[10];
    if (Math.abs(h - u) < 0.01 && Math.abs(f - x) < 0.01 && Math.abs(g - p) < 0.01) {
      if (Math.abs(h + u) < 0.1 && Math.abs(f + x) < 0.1 && Math.abs(g + p) < 0.1 && Math.abs(c + m + d - 3) < 0.1)
        return this.set(1, 0, 0, 0), this;
      t = Math.PI;
      const _ = (c + 1) / 2, T = (m + 1) / 2, C = (d + 1) / 2, L = (h + u) / 4, w = (f + x) / 4, V = (g + p) / 4;
      return _ > T && _ > C ? _ < 0.01 ? (n = 0, r = 0.707106781, s = 0.707106781) : (n = Math.sqrt(_), r = L / n, s = w / n) : T > C ? T < 0.01 ? (n = 0.707106781, r = 0, s = 0.707106781) : (r = Math.sqrt(T), n = L / r, s = V / r) : C < 0.01 ? (n = 0.707106781, r = 0.707106781, s = 0) : (s = Math.sqrt(C), n = w / s, r = V / s), this.set(n, r, s, t), this;
    }
    let A = Math.sqrt((p - g) * (p - g) + (f - x) * (f - x) + (u - h) * (u - h));
    return Math.abs(A) < 1e-3 && (A = 1), this.x = (p - g) / A, this.y = (f - x) / A, this.z = (u - h) / A, this.w = Math.acos((c + m + d - 1) / 2), this;
  }
  min(e) {
    return this.x = Math.min(this.x, e.x), this.y = Math.min(this.y, e.y), this.z = Math.min(this.z, e.z), this.w = Math.min(this.w, e.w), this;
  }
  max(e) {
    return this.x = Math.max(this.x, e.x), this.y = Math.max(this.y, e.y), this.z = Math.max(this.z, e.z), this.w = Math.max(this.w, e.w), this;
  }
  clamp(e, t) {
    return this.x = Math.max(e.x, Math.min(t.x, this.x)), this.y = Math.max(e.y, Math.min(t.y, this.y)), this.z = Math.max(e.z, Math.min(t.z, this.z)), this.w = Math.max(e.w, Math.min(t.w, this.w)), this;
  }
  clampScalar(e, t) {
    return this.x = Math.max(e, Math.min(t, this.x)), this.y = Math.max(e, Math.min(t, this.y)), this.z = Math.max(e, Math.min(t, this.z)), this.w = Math.max(e, Math.min(t, this.w)), this;
  }
  clampLength(e, t) {
    const n = this.length();
    return this.divideScalar(n || 1).multiplyScalar(Math.max(e, Math.min(t, n)));
  }
  floor() {
    return this.x = Math.floor(this.x), this.y = Math.floor(this.y), this.z = Math.floor(this.z), this.w = Math.floor(this.w), this;
  }
  ceil() {
    return this.x = Math.ceil(this.x), this.y = Math.ceil(this.y), this.z = Math.ceil(this.z), this.w = Math.ceil(this.w), this;
  }
  round() {
    return this.x = Math.round(this.x), this.y = Math.round(this.y), this.z = Math.round(this.z), this.w = Math.round(this.w), this;
  }
  roundToZero() {
    return this.x = this.x < 0 ? Math.ceil(this.x) : Math.floor(this.x), this.y = this.y < 0 ? Math.ceil(this.y) : Math.floor(this.y), this.z = this.z < 0 ? Math.ceil(this.z) : Math.floor(this.z), this.w = this.w < 0 ? Math.ceil(this.w) : Math.floor(this.w), this;
  }
  negate() {
    return this.x = -this.x, this.y = -this.y, this.z = -this.z, this.w = -this.w, this;
  }
  dot(e) {
    return this.x * e.x + this.y * e.y + this.z * e.z + this.w * e.w;
  }
  lengthSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w;
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z + this.w * this.w);
  }
  manhattanLength() {
    return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z) + Math.abs(this.w);
  }
  normalize() {
    return this.divideScalar(this.length() || 1);
  }
  setLength(e) {
    return this.normalize().multiplyScalar(e);
  }
  lerp(e, t) {
    return this.x += (e.x - this.x) * t, this.y += (e.y - this.y) * t, this.z += (e.z - this.z) * t, this.w += (e.w - this.w) * t, this;
  }
  lerpVectors(e, t, n) {
    return this.x = e.x + (t.x - e.x) * n, this.y = e.y + (t.y - e.y) * n, this.z = e.z + (t.z - e.z) * n, this.w = e.w + (t.w - e.w) * n, this;
  }
  equals(e) {
    return e.x === this.x && e.y === this.y && e.z === this.z && e.w === this.w;
  }
  fromArray(e, t = 0) {
    return this.x = e[t], this.y = e[t + 1], this.z = e[t + 2], this.w = e[t + 3], this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this.x, e[t + 1] = this.y, e[t + 2] = this.z, e[t + 3] = this.w, e;
  }
  fromBufferAttribute(e, t) {
    return this.x = e.getX(t), this.y = e.getY(t), this.z = e.getZ(t), this.w = e.getW(t), this;
  }
  random() {
    return this.x = Math.random(), this.y = Math.random(), this.z = Math.random(), this.w = Math.random(), this;
  }
  *[Symbol.iterator]() {
    yield this.x, yield this.y, yield this.z, yield this.w;
  }
}
class mc extends Rn {
  constructor(e = 1, t = 1, n = {}) {
    super(), this.isRenderTarget = !0, this.width = e, this.height = t, this.depth = 1, this.scissor = new ot(0, 0, e, t), this.scissorTest = !1, this.viewport = new ot(0, 0, e, t);
    const r = { width: e, height: t, depth: 1 };
    n.encoding !== void 0 && (di("THREE.WebGLRenderTarget: option.encoding has been replaced by option.colorSpace."), n.colorSpace = n.encoding === Sn ? Oe : En), this.texture = new St(r, n.mapping, n.wrapS, n.wrapT, n.magFilter, n.minFilter, n.format, n.type, n.anisotropy, n.colorSpace), this.texture.isRenderTargetTexture = !0, this.texture.flipY = !1, this.texture.generateMipmaps = n.generateMipmaps !== void 0 ? n.generateMipmaps : !1, this.texture.internalFormat = n.internalFormat !== void 0 ? n.internalFormat : null, this.texture.minFilter = n.minFilter !== void 0 ? n.minFilter : Tt, this.depthBuffer = n.depthBuffer !== void 0 ? n.depthBuffer : !0, this.stencilBuffer = n.stencilBuffer !== void 0 ? n.stencilBuffer : !1, this.depthTexture = n.depthTexture !== void 0 ? n.depthTexture : null, this.samples = n.samples !== void 0 ? n.samples : 0;
  }
  setSize(e, t, n = 1) {
    (this.width !== e || this.height !== t || this.depth !== n) && (this.width = e, this.height = t, this.depth = n, this.texture.image.width = e, this.texture.image.height = t, this.texture.image.depth = n, this.dispose()), this.viewport.set(0, 0, e, t), this.scissor.set(0, 0, e, t);
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    this.width = e.width, this.height = e.height, this.depth = e.depth, this.scissor.copy(e.scissor), this.scissorTest = e.scissorTest, this.viewport.copy(e.viewport), this.texture = e.texture.clone(), this.texture.isRenderTargetTexture = !0;
    const t = Object.assign({}, e.texture.image);
    return this.texture.source = new mo(t), this.depthBuffer = e.depthBuffer, this.stencilBuffer = e.stencilBuffer, e.depthTexture !== null && (this.depthTexture = e.depthTexture.clone()), this.samples = e.samples, this;
  }
  dispose() {
    this.dispatchEvent({ type: "dispose" });
  }
}
class Tn extends mc {
  constructor(e = 1, t = 1, n = {}) {
    super(e, t, n), this.isWebGLRenderTarget = !0;
  }
}
class go extends St {
  constructor(e = null, t = 1, n = 1, r = 1) {
    super(null), this.isDataArrayTexture = !0, this.image = { data: e, width: t, height: n, depth: r }, this.magFilter = mt, this.minFilter = mt, this.wrapR = Lt, this.generateMipmaps = !1, this.flipY = !1, this.unpackAlignment = 1;
  }
}
class gc extends St {
  constructor(e = null, t = 1, n = 1, r = 1) {
    super(null), this.isData3DTexture = !0, this.image = { data: e, width: t, height: n, depth: r }, this.magFilter = mt, this.minFilter = mt, this.wrapR = Lt, this.generateMipmaps = !1, this.flipY = !1, this.unpackAlignment = 1;
  }
}
class An {
  constructor(e = 0, t = 0, n = 0, r = 1) {
    this.isQuaternion = !0, this._x = e, this._y = t, this._z = n, this._w = r;
  }
  static slerpFlat(e, t, n, r, s, o, a) {
    let l = n[r + 0], c = n[r + 1], h = n[r + 2], f = n[r + 3];
    const u = s[o + 0], m = s[o + 1], g = s[o + 2], x = s[o + 3];
    if (a === 0) {
      e[t + 0] = l, e[t + 1] = c, e[t + 2] = h, e[t + 3] = f;
      return;
    }
    if (a === 1) {
      e[t + 0] = u, e[t + 1] = m, e[t + 2] = g, e[t + 3] = x;
      return;
    }
    if (f !== x || l !== u || c !== m || h !== g) {
      let p = 1 - a;
      const d = l * u + c * m + h * g + f * x, A = d >= 0 ? 1 : -1, _ = 1 - d * d;
      if (_ > Number.EPSILON) {
        const C = Math.sqrt(_), L = Math.atan2(C, d * A);
        p = Math.sin(p * L) / C, a = Math.sin(a * L) / C;
      }
      const T = a * A;
      if (l = l * p + u * T, c = c * p + m * T, h = h * p + g * T, f = f * p + x * T, p === 1 - a) {
        const C = 1 / Math.sqrt(l * l + c * c + h * h + f * f);
        l *= C, c *= C, h *= C, f *= C;
      }
    }
    e[t] = l, e[t + 1] = c, e[t + 2] = h, e[t + 3] = f;
  }
  static multiplyQuaternionsFlat(e, t, n, r, s, o) {
    const a = n[r], l = n[r + 1], c = n[r + 2], h = n[r + 3], f = s[o], u = s[o + 1], m = s[o + 2], g = s[o + 3];
    return e[t] = a * g + h * f + l * m - c * u, e[t + 1] = l * g + h * u + c * f - a * m, e[t + 2] = c * g + h * m + a * u - l * f, e[t + 3] = h * g - a * f - l * u - c * m, e;
  }
  get x() {
    return this._x;
  }
  set x(e) {
    this._x = e, this._onChangeCallback();
  }
  get y() {
    return this._y;
  }
  set y(e) {
    this._y = e, this._onChangeCallback();
  }
  get z() {
    return this._z;
  }
  set z(e) {
    this._z = e, this._onChangeCallback();
  }
  get w() {
    return this._w;
  }
  set w(e) {
    this._w = e, this._onChangeCallback();
  }
  set(e, t, n, r) {
    return this._x = e, this._y = t, this._z = n, this._w = r, this._onChangeCallback(), this;
  }
  clone() {
    return new this.constructor(this._x, this._y, this._z, this._w);
  }
  copy(e) {
    return this._x = e.x, this._y = e.y, this._z = e.z, this._w = e.w, this._onChangeCallback(), this;
  }
  setFromEuler(e, t) {
    const n = e._x, r = e._y, s = e._z, o = e._order, a = Math.cos, l = Math.sin, c = a(n / 2), h = a(r / 2), f = a(s / 2), u = l(n / 2), m = l(r / 2), g = l(s / 2);
    switch (o) {
      case "XYZ":
        this._x = u * h * f + c * m * g, this._y = c * m * f - u * h * g, this._z = c * h * g + u * m * f, this._w = c * h * f - u * m * g;
        break;
      case "YXZ":
        this._x = u * h * f + c * m * g, this._y = c * m * f - u * h * g, this._z = c * h * g - u * m * f, this._w = c * h * f + u * m * g;
        break;
      case "ZXY":
        this._x = u * h * f - c * m * g, this._y = c * m * f + u * h * g, this._z = c * h * g + u * m * f, this._w = c * h * f - u * m * g;
        break;
      case "ZYX":
        this._x = u * h * f - c * m * g, this._y = c * m * f + u * h * g, this._z = c * h * g - u * m * f, this._w = c * h * f + u * m * g;
        break;
      case "YZX":
        this._x = u * h * f + c * m * g, this._y = c * m * f + u * h * g, this._z = c * h * g - u * m * f, this._w = c * h * f - u * m * g;
        break;
      case "XZY":
        this._x = u * h * f - c * m * g, this._y = c * m * f - u * h * g, this._z = c * h * g + u * m * f, this._w = c * h * f + u * m * g;
        break;
      default:
        console.warn("THREE.Quaternion: .setFromEuler() encountered an unknown order: " + o);
    }
    return t !== !1 && this._onChangeCallback(), this;
  }
  setFromAxisAngle(e, t) {
    const n = t / 2, r = Math.sin(n);
    return this._x = e.x * r, this._y = e.y * r, this._z = e.z * r, this._w = Math.cos(n), this._onChangeCallback(), this;
  }
  setFromRotationMatrix(e) {
    const t = e.elements, n = t[0], r = t[4], s = t[8], o = t[1], a = t[5], l = t[9], c = t[2], h = t[6], f = t[10], u = n + a + f;
    if (u > 0) {
      const m = 0.5 / Math.sqrt(u + 1);
      this._w = 0.25 / m, this._x = (h - l) * m, this._y = (s - c) * m, this._z = (o - r) * m;
    } else if (n > a && n > f) {
      const m = 2 * Math.sqrt(1 + n - a - f);
      this._w = (h - l) / m, this._x = 0.25 * m, this._y = (r + o) / m, this._z = (s + c) / m;
    } else if (a > f) {
      const m = 2 * Math.sqrt(1 + a - n - f);
      this._w = (s - c) / m, this._x = (r + o) / m, this._y = 0.25 * m, this._z = (l + h) / m;
    } else {
      const m = 2 * Math.sqrt(1 + f - n - a);
      this._w = (o - r) / m, this._x = (s + c) / m, this._y = (l + h) / m, this._z = 0.25 * m;
    }
    return this._onChangeCallback(), this;
  }
  setFromUnitVectors(e, t) {
    let n = e.dot(t) + 1;
    return n < Number.EPSILON ? (n = 0, Math.abs(e.x) > Math.abs(e.z) ? (this._x = -e.y, this._y = e.x, this._z = 0, this._w = n) : (this._x = 0, this._y = -e.z, this._z = e.y, this._w = n)) : (this._x = e.y * t.z - e.z * t.y, this._y = e.z * t.x - e.x * t.z, this._z = e.x * t.y - e.y * t.x, this._w = n), this.normalize();
  }
  angleTo(e) {
    return 2 * Math.acos(Math.abs(it(this.dot(e), -1, 1)));
  }
  rotateTowards(e, t) {
    const n = this.angleTo(e);
    if (n === 0)
      return this;
    const r = Math.min(1, t / n);
    return this.slerp(e, r), this;
  }
  identity() {
    return this.set(0, 0, 0, 1);
  }
  invert() {
    return this.conjugate();
  }
  conjugate() {
    return this._x *= -1, this._y *= -1, this._z *= -1, this._onChangeCallback(), this;
  }
  dot(e) {
    return this._x * e._x + this._y * e._y + this._z * e._z + this._w * e._w;
  }
  lengthSq() {
    return this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w;
  }
  length() {
    return Math.sqrt(this._x * this._x + this._y * this._y + this._z * this._z + this._w * this._w);
  }
  normalize() {
    let e = this.length();
    return e === 0 ? (this._x = 0, this._y = 0, this._z = 0, this._w = 1) : (e = 1 / e, this._x = this._x * e, this._y = this._y * e, this._z = this._z * e, this._w = this._w * e), this._onChangeCallback(), this;
  }
  multiply(e) {
    return this.multiplyQuaternions(this, e);
  }
  premultiply(e) {
    return this.multiplyQuaternions(e, this);
  }
  multiplyQuaternions(e, t) {
    const n = e._x, r = e._y, s = e._z, o = e._w, a = t._x, l = t._y, c = t._z, h = t._w;
    return this._x = n * h + o * a + r * c - s * l, this._y = r * h + o * l + s * a - n * c, this._z = s * h + o * c + n * l - r * a, this._w = o * h - n * a - r * l - s * c, this._onChangeCallback(), this;
  }
  slerp(e, t) {
    if (t === 0)
      return this;
    if (t === 1)
      return this.copy(e);
    const n = this._x, r = this._y, s = this._z, o = this._w;
    let a = o * e._w + n * e._x + r * e._y + s * e._z;
    if (a < 0 ? (this._w = -e._w, this._x = -e._x, this._y = -e._y, this._z = -e._z, a = -a) : this.copy(e), a >= 1)
      return this._w = o, this._x = n, this._y = r, this._z = s, this;
    const l = 1 - a * a;
    if (l <= Number.EPSILON) {
      const m = 1 - t;
      return this._w = m * o + t * this._w, this._x = m * n + t * this._x, this._y = m * r + t * this._y, this._z = m * s + t * this._z, this.normalize(), this._onChangeCallback(), this;
    }
    const c = Math.sqrt(l), h = Math.atan2(c, a), f = Math.sin((1 - t) * h) / c, u = Math.sin(t * h) / c;
    return this._w = o * f + this._w * u, this._x = n * f + this._x * u, this._y = r * f + this._y * u, this._z = s * f + this._z * u, this._onChangeCallback(), this;
  }
  slerpQuaternions(e, t, n) {
    return this.copy(e).slerp(t, n);
  }
  random() {
    const e = Math.random(), t = Math.sqrt(1 - e), n = Math.sqrt(e), r = 2 * Math.PI * Math.random(), s = 2 * Math.PI * Math.random();
    return this.set(
      t * Math.cos(r),
      n * Math.sin(s),
      n * Math.cos(s),
      t * Math.sin(r)
    );
  }
  equals(e) {
    return e._x === this._x && e._y === this._y && e._z === this._z && e._w === this._w;
  }
  fromArray(e, t = 0) {
    return this._x = e[t], this._y = e[t + 1], this._z = e[t + 2], this._w = e[t + 3], this._onChangeCallback(), this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this._x, e[t + 1] = this._y, e[t + 2] = this._z, e[t + 3] = this._w, e;
  }
  fromBufferAttribute(e, t) {
    return this._x = e.getX(t), this._y = e.getY(t), this._z = e.getZ(t), this._w = e.getW(t), this;
  }
  toJSON() {
    return this.toArray();
  }
  _onChange(e) {
    return this._onChangeCallback = e, this;
  }
  _onChangeCallback() {
  }
  *[Symbol.iterator]() {
    yield this._x, yield this._y, yield this._z, yield this._w;
  }
}
class U {
  constructor(e = 0, t = 0, n = 0) {
    U.prototype.isVector3 = !0, this.x = e, this.y = t, this.z = n;
  }
  set(e, t, n) {
    return n === void 0 && (n = this.z), this.x = e, this.y = t, this.z = n, this;
  }
  setScalar(e) {
    return this.x = e, this.y = e, this.z = e, this;
  }
  setX(e) {
    return this.x = e, this;
  }
  setY(e) {
    return this.y = e, this;
  }
  setZ(e) {
    return this.z = e, this;
  }
  setComponent(e, t) {
    switch (e) {
      case 0:
        this.x = t;
        break;
      case 1:
        this.y = t;
        break;
      case 2:
        this.z = t;
        break;
      default:
        throw new Error("index is out of range: " + e);
    }
    return this;
  }
  getComponent(e) {
    switch (e) {
      case 0:
        return this.x;
      case 1:
        return this.y;
      case 2:
        return this.z;
      default:
        throw new Error("index is out of range: " + e);
    }
  }
  clone() {
    return new this.constructor(this.x, this.y, this.z);
  }
  copy(e) {
    return this.x = e.x, this.y = e.y, this.z = e.z, this;
  }
  add(e) {
    return this.x += e.x, this.y += e.y, this.z += e.z, this;
  }
  addScalar(e) {
    return this.x += e, this.y += e, this.z += e, this;
  }
  addVectors(e, t) {
    return this.x = e.x + t.x, this.y = e.y + t.y, this.z = e.z + t.z, this;
  }
  addScaledVector(e, t) {
    return this.x += e.x * t, this.y += e.y * t, this.z += e.z * t, this;
  }
  sub(e) {
    return this.x -= e.x, this.y -= e.y, this.z -= e.z, this;
  }
  subScalar(e) {
    return this.x -= e, this.y -= e, this.z -= e, this;
  }
  subVectors(e, t) {
    return this.x = e.x - t.x, this.y = e.y - t.y, this.z = e.z - t.z, this;
  }
  multiply(e) {
    return this.x *= e.x, this.y *= e.y, this.z *= e.z, this;
  }
  multiplyScalar(e) {
    return this.x *= e, this.y *= e, this.z *= e, this;
  }
  multiplyVectors(e, t) {
    return this.x = e.x * t.x, this.y = e.y * t.y, this.z = e.z * t.z, this;
  }
  applyEuler(e) {
    return this.applyQuaternion(ia.setFromEuler(e));
  }
  applyAxisAngle(e, t) {
    return this.applyQuaternion(ia.setFromAxisAngle(e, t));
  }
  applyMatrix3(e) {
    const t = this.x, n = this.y, r = this.z, s = e.elements;
    return this.x = s[0] * t + s[3] * n + s[6] * r, this.y = s[1] * t + s[4] * n + s[7] * r, this.z = s[2] * t + s[5] * n + s[8] * r, this;
  }
  applyNormalMatrix(e) {
    return this.applyMatrix3(e).normalize();
  }
  applyMatrix4(e) {
    const t = this.x, n = this.y, r = this.z, s = e.elements, o = 1 / (s[3] * t + s[7] * n + s[11] * r + s[15]);
    return this.x = (s[0] * t + s[4] * n + s[8] * r + s[12]) * o, this.y = (s[1] * t + s[5] * n + s[9] * r + s[13]) * o, this.z = (s[2] * t + s[6] * n + s[10] * r + s[14]) * o, this;
  }
  applyQuaternion(e) {
    const t = this.x, n = this.y, r = this.z, s = e.x, o = e.y, a = e.z, l = e.w, c = l * t + o * r - a * n, h = l * n + a * t - s * r, f = l * r + s * n - o * t, u = -s * t - o * n - a * r;
    return this.x = c * l + u * -s + h * -a - f * -o, this.y = h * l + u * -o + f * -s - c * -a, this.z = f * l + u * -a + c * -o - h * -s, this;
  }
  project(e) {
    return this.applyMatrix4(e.matrixWorldInverse).applyMatrix4(e.projectionMatrix);
  }
  unproject(e) {
    return this.applyMatrix4(e.projectionMatrixInverse).applyMatrix4(e.matrixWorld);
  }
  transformDirection(e) {
    const t = this.x, n = this.y, r = this.z, s = e.elements;
    return this.x = s[0] * t + s[4] * n + s[8] * r, this.y = s[1] * t + s[5] * n + s[9] * r, this.z = s[2] * t + s[6] * n + s[10] * r, this.normalize();
  }
  divide(e) {
    return this.x /= e.x, this.y /= e.y, this.z /= e.z, this;
  }
  divideScalar(e) {
    return this.multiplyScalar(1 / e);
  }
  min(e) {
    return this.x = Math.min(this.x, e.x), this.y = Math.min(this.y, e.y), this.z = Math.min(this.z, e.z), this;
  }
  max(e) {
    return this.x = Math.max(this.x, e.x), this.y = Math.max(this.y, e.y), this.z = Math.max(this.z, e.z), this;
  }
  clamp(e, t) {
    return this.x = Math.max(e.x, Math.min(t.x, this.x)), this.y = Math.max(e.y, Math.min(t.y, this.y)), this.z = Math.max(e.z, Math.min(t.z, this.z)), this;
  }
  clampScalar(e, t) {
    return this.x = Math.max(e, Math.min(t, this.x)), this.y = Math.max(e, Math.min(t, this.y)), this.z = Math.max(e, Math.min(t, this.z)), this;
  }
  clampLength(e, t) {
    const n = this.length();
    return this.divideScalar(n || 1).multiplyScalar(Math.max(e, Math.min(t, n)));
  }
  floor() {
    return this.x = Math.floor(this.x), this.y = Math.floor(this.y), this.z = Math.floor(this.z), this;
  }
  ceil() {
    return this.x = Math.ceil(this.x), this.y = Math.ceil(this.y), this.z = Math.ceil(this.z), this;
  }
  round() {
    return this.x = Math.round(this.x), this.y = Math.round(this.y), this.z = Math.round(this.z), this;
  }
  roundToZero() {
    return this.x = this.x < 0 ? Math.ceil(this.x) : Math.floor(this.x), this.y = this.y < 0 ? Math.ceil(this.y) : Math.floor(this.y), this.z = this.z < 0 ? Math.ceil(this.z) : Math.floor(this.z), this;
  }
  negate() {
    return this.x = -this.x, this.y = -this.y, this.z = -this.z, this;
  }
  dot(e) {
    return this.x * e.x + this.y * e.y + this.z * e.z;
  }
  // TODO lengthSquared?
  lengthSq() {
    return this.x * this.x + this.y * this.y + this.z * this.z;
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
  }
  manhattanLength() {
    return Math.abs(this.x) + Math.abs(this.y) + Math.abs(this.z);
  }
  normalize() {
    return this.divideScalar(this.length() || 1);
  }
  setLength(e) {
    return this.normalize().multiplyScalar(e);
  }
  lerp(e, t) {
    return this.x += (e.x - this.x) * t, this.y += (e.y - this.y) * t, this.z += (e.z - this.z) * t, this;
  }
  lerpVectors(e, t, n) {
    return this.x = e.x + (t.x - e.x) * n, this.y = e.y + (t.y - e.y) * n, this.z = e.z + (t.z - e.z) * n, this;
  }
  cross(e) {
    return this.crossVectors(this, e);
  }
  crossVectors(e, t) {
    const n = e.x, r = e.y, s = e.z, o = t.x, a = t.y, l = t.z;
    return this.x = r * l - s * a, this.y = s * o - n * l, this.z = n * a - r * o, this;
  }
  projectOnVector(e) {
    const t = e.lengthSq();
    if (t === 0)
      return this.set(0, 0, 0);
    const n = e.dot(this) / t;
    return this.copy(e).multiplyScalar(n);
  }
  projectOnPlane(e) {
    return yr.copy(this).projectOnVector(e), this.sub(yr);
  }
  reflect(e) {
    return this.sub(yr.copy(e).multiplyScalar(2 * this.dot(e)));
  }
  angleTo(e) {
    const t = Math.sqrt(this.lengthSq() * e.lengthSq());
    if (t === 0)
      return Math.PI / 2;
    const n = this.dot(e) / t;
    return Math.acos(it(n, -1, 1));
  }
  distanceTo(e) {
    return Math.sqrt(this.distanceToSquared(e));
  }
  distanceToSquared(e) {
    const t = this.x - e.x, n = this.y - e.y, r = this.z - e.z;
    return t * t + n * n + r * r;
  }
  manhattanDistanceTo(e) {
    return Math.abs(this.x - e.x) + Math.abs(this.y - e.y) + Math.abs(this.z - e.z);
  }
  setFromSpherical(e) {
    return this.setFromSphericalCoords(e.radius, e.phi, e.theta);
  }
  setFromSphericalCoords(e, t, n) {
    const r = Math.sin(t) * e;
    return this.x = r * Math.sin(n), this.y = Math.cos(t) * e, this.z = r * Math.cos(n), this;
  }
  setFromCylindrical(e) {
    return this.setFromCylindricalCoords(e.radius, e.theta, e.y);
  }
  setFromCylindricalCoords(e, t, n) {
    return this.x = e * Math.sin(t), this.y = n, this.z = e * Math.cos(t), this;
  }
  setFromMatrixPosition(e) {
    const t = e.elements;
    return this.x = t[12], this.y = t[13], this.z = t[14], this;
  }
  setFromMatrixScale(e) {
    const t = this.setFromMatrixColumn(e, 0).length(), n = this.setFromMatrixColumn(e, 1).length(), r = this.setFromMatrixColumn(e, 2).length();
    return this.x = t, this.y = n, this.z = r, this;
  }
  setFromMatrixColumn(e, t) {
    return this.fromArray(e.elements, t * 4);
  }
  setFromMatrix3Column(e, t) {
    return this.fromArray(e.elements, t * 3);
  }
  setFromEuler(e) {
    return this.x = e._x, this.y = e._y, this.z = e._z, this;
  }
  setFromColor(e) {
    return this.x = e.r, this.y = e.g, this.z = e.b, this;
  }
  equals(e) {
    return e.x === this.x && e.y === this.y && e.z === this.z;
  }
  fromArray(e, t = 0) {
    return this.x = e[t], this.y = e[t + 1], this.z = e[t + 2], this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this.x, e[t + 1] = this.y, e[t + 2] = this.z, e;
  }
  fromBufferAttribute(e, t) {
    return this.x = e.getX(t), this.y = e.getY(t), this.z = e.getZ(t), this;
  }
  random() {
    return this.x = Math.random(), this.y = Math.random(), this.z = Math.random(), this;
  }
  randomDirection() {
    const e = (Math.random() - 0.5) * 2, t = Math.random() * Math.PI * 2, n = Math.sqrt(1 - e ** 2);
    return this.x = n * Math.cos(t), this.y = n * Math.sin(t), this.z = e, this;
  }
  *[Symbol.iterator]() {
    yield this.x, yield this.y, yield this.z;
  }
}
const yr = /* @__PURE__ */ new U(), ia = /* @__PURE__ */ new An();
class Ti {
  constructor(e = new U(1 / 0, 1 / 0, 1 / 0), t = new U(-1 / 0, -1 / 0, -1 / 0)) {
    this.isBox3 = !0, this.min = e, this.max = t;
  }
  set(e, t) {
    return this.min.copy(e), this.max.copy(t), this;
  }
  setFromArray(e) {
    this.makeEmpty();
    for (let t = 0, n = e.length; t < n; t += 3)
      this.expandByPoint(Gt.fromArray(e, t));
    return this;
  }
  setFromBufferAttribute(e) {
    this.makeEmpty();
    for (let t = 0, n = e.count; t < n; t++)
      this.expandByPoint(Gt.fromBufferAttribute(e, t));
    return this;
  }
  setFromPoints(e) {
    this.makeEmpty();
    for (let t = 0, n = e.length; t < n; t++)
      this.expandByPoint(e[t]);
    return this;
  }
  setFromCenterAndSize(e, t) {
    const n = Gt.copy(t).multiplyScalar(0.5);
    return this.min.copy(e).sub(n), this.max.copy(e).add(n), this;
  }
  setFromObject(e, t = !1) {
    return this.makeEmpty(), this.expandByObject(e, t);
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    return this.min.copy(e.min), this.max.copy(e.max), this;
  }
  makeEmpty() {
    return this.min.x = this.min.y = this.min.z = 1 / 0, this.max.x = this.max.y = this.max.z = -1 / 0, this;
  }
  isEmpty() {
    return this.max.x < this.min.x || this.max.y < this.min.y || this.max.z < this.min.z;
  }
  getCenter(e) {
    return this.isEmpty() ? e.set(0, 0, 0) : e.addVectors(this.min, this.max).multiplyScalar(0.5);
  }
  getSize(e) {
    return this.isEmpty() ? e.set(0, 0, 0) : e.subVectors(this.max, this.min);
  }
  expandByPoint(e) {
    return this.min.min(e), this.max.max(e), this;
  }
  expandByVector(e) {
    return this.min.sub(e), this.max.add(e), this;
  }
  expandByScalar(e) {
    return this.min.addScalar(-e), this.max.addScalar(e), this;
  }
  expandByObject(e, t = !1) {
    if (e.updateWorldMatrix(!1, !1), e.boundingBox !== void 0)
      e.boundingBox === null && e.computeBoundingBox(), Dn.copy(e.boundingBox), Dn.applyMatrix4(e.matrixWorld), this.union(Dn);
    else {
      const r = e.geometry;
      if (r !== void 0)
        if (t && r.attributes !== void 0 && r.attributes.position !== void 0) {
          const s = r.attributes.position;
          for (let o = 0, a = s.count; o < a; o++)
            Gt.fromBufferAttribute(s, o).applyMatrix4(e.matrixWorld), this.expandByPoint(Gt);
        } else
          r.boundingBox === null && r.computeBoundingBox(), Dn.copy(r.boundingBox), Dn.applyMatrix4(e.matrixWorld), this.union(Dn);
    }
    const n = e.children;
    for (let r = 0, s = n.length; r < s; r++)
      this.expandByObject(n[r], t);
    return this;
  }
  containsPoint(e) {
    return !(e.x < this.min.x || e.x > this.max.x || e.y < this.min.y || e.y > this.max.y || e.z < this.min.z || e.z > this.max.z);
  }
  containsBox(e) {
    return this.min.x <= e.min.x && e.max.x <= this.max.x && this.min.y <= e.min.y && e.max.y <= this.max.y && this.min.z <= e.min.z && e.max.z <= this.max.z;
  }
  getParameter(e, t) {
    return t.set(
      (e.x - this.min.x) / (this.max.x - this.min.x),
      (e.y - this.min.y) / (this.max.y - this.min.y),
      (e.z - this.min.z) / (this.max.z - this.min.z)
    );
  }
  intersectsBox(e) {
    return !(e.max.x < this.min.x || e.min.x > this.max.x || e.max.y < this.min.y || e.min.y > this.max.y || e.max.z < this.min.z || e.min.z > this.max.z);
  }
  intersectsSphere(e) {
    return this.clampPoint(e.center, Gt), Gt.distanceToSquared(e.center) <= e.radius * e.radius;
  }
  intersectsPlane(e) {
    let t, n;
    return e.normal.x > 0 ? (t = e.normal.x * this.min.x, n = e.normal.x * this.max.x) : (t = e.normal.x * this.max.x, n = e.normal.x * this.min.x), e.normal.y > 0 ? (t += e.normal.y * this.min.y, n += e.normal.y * this.max.y) : (t += e.normal.y * this.max.y, n += e.normal.y * this.min.y), e.normal.z > 0 ? (t += e.normal.z * this.min.z, n += e.normal.z * this.max.z) : (t += e.normal.z * this.max.z, n += e.normal.z * this.min.z), t <= -e.constant && n >= -e.constant;
  }
  intersectsTriangle(e) {
    if (this.isEmpty())
      return !1;
    this.getCenter(ai), Ci.subVectors(this.max, ai), In.subVectors(e.a, ai), Nn.subVectors(e.b, ai), On.subVectors(e.c, ai), Kt.subVectors(Nn, In), Jt.subVectors(On, Nn), fn.subVectors(In, On);
    let t = [
      0,
      -Kt.z,
      Kt.y,
      0,
      -Jt.z,
      Jt.y,
      0,
      -fn.z,
      fn.y,
      Kt.z,
      0,
      -Kt.x,
      Jt.z,
      0,
      -Jt.x,
      fn.z,
      0,
      -fn.x,
      -Kt.y,
      Kt.x,
      0,
      -Jt.y,
      Jt.x,
      0,
      -fn.y,
      fn.x,
      0
    ];
    return !Tr(t, In, Nn, On, Ci) || (t = [1, 0, 0, 0, 1, 0, 0, 0, 1], !Tr(t, In, Nn, On, Ci)) ? !1 : (Pi.crossVectors(Kt, Jt), t = [Pi.x, Pi.y, Pi.z], Tr(t, In, Nn, On, Ci));
  }
  clampPoint(e, t) {
    return t.copy(e).clamp(this.min, this.max);
  }
  distanceToPoint(e) {
    return this.clampPoint(e, Gt).distanceTo(e);
  }
  getBoundingSphere(e) {
    return this.isEmpty() ? e.makeEmpty() : (this.getCenter(e.center), e.radius = this.getSize(Gt).length() * 0.5), e;
  }
  intersect(e) {
    return this.min.max(e.min), this.max.min(e.max), this.isEmpty() && this.makeEmpty(), this;
  }
  union(e) {
    return this.min.min(e.min), this.max.max(e.max), this;
  }
  applyMatrix4(e) {
    return this.isEmpty() ? this : (Ht[0].set(this.min.x, this.min.y, this.min.z).applyMatrix4(e), Ht[1].set(this.min.x, this.min.y, this.max.z).applyMatrix4(e), Ht[2].set(this.min.x, this.max.y, this.min.z).applyMatrix4(e), Ht[3].set(this.min.x, this.max.y, this.max.z).applyMatrix4(e), Ht[4].set(this.max.x, this.min.y, this.min.z).applyMatrix4(e), Ht[5].set(this.max.x, this.min.y, this.max.z).applyMatrix4(e), Ht[6].set(this.max.x, this.max.y, this.min.z).applyMatrix4(e), Ht[7].set(this.max.x, this.max.y, this.max.z).applyMatrix4(e), this.setFromPoints(Ht), this);
  }
  translate(e) {
    return this.min.add(e), this.max.add(e), this;
  }
  equals(e) {
    return e.min.equals(this.min) && e.max.equals(this.max);
  }
}
const Ht = [
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U(),
  /* @__PURE__ */ new U()
], Gt = /* @__PURE__ */ new U(), Dn = /* @__PURE__ */ new Ti(), In = /* @__PURE__ */ new U(), Nn = /* @__PURE__ */ new U(), On = /* @__PURE__ */ new U(), Kt = /* @__PURE__ */ new U(), Jt = /* @__PURE__ */ new U(), fn = /* @__PURE__ */ new U(), ai = /* @__PURE__ */ new U(), Ci = /* @__PURE__ */ new U(), Pi = /* @__PURE__ */ new U(), dn = /* @__PURE__ */ new U();
function Tr(i, e, t, n, r) {
  for (let s = 0, o = i.length - 3; s <= o; s += 3) {
    dn.fromArray(i, s);
    const a = r.x * Math.abs(dn.x) + r.y * Math.abs(dn.y) + r.z * Math.abs(dn.z), l = e.dot(dn), c = t.dot(dn), h = n.dot(dn);
    if (Math.max(-Math.max(l, c, h), Math.min(l, c, h)) > a)
      return !1;
  }
  return !0;
}
const _c = /* @__PURE__ */ new Ti(), oi = /* @__PURE__ */ new U(), Ar = /* @__PURE__ */ new U();
class ls {
  constructor(e = new U(), t = -1) {
    this.center = e, this.radius = t;
  }
  set(e, t) {
    return this.center.copy(e), this.radius = t, this;
  }
  setFromPoints(e, t) {
    const n = this.center;
    t !== void 0 ? n.copy(t) : _c.setFromPoints(e).getCenter(n);
    let r = 0;
    for (let s = 0, o = e.length; s < o; s++)
      r = Math.max(r, n.distanceToSquared(e[s]));
    return this.radius = Math.sqrt(r), this;
  }
  copy(e) {
    return this.center.copy(e.center), this.radius = e.radius, this;
  }
  isEmpty() {
    return this.radius < 0;
  }
  makeEmpty() {
    return this.center.set(0, 0, 0), this.radius = -1, this;
  }
  containsPoint(e) {
    return e.distanceToSquared(this.center) <= this.radius * this.radius;
  }
  distanceToPoint(e) {
    return e.distanceTo(this.center) - this.radius;
  }
  intersectsSphere(e) {
    const t = this.radius + e.radius;
    return e.center.distanceToSquared(this.center) <= t * t;
  }
  intersectsBox(e) {
    return e.intersectsSphere(this);
  }
  intersectsPlane(e) {
    return Math.abs(e.distanceToPoint(this.center)) <= this.radius;
  }
  clampPoint(e, t) {
    const n = this.center.distanceToSquared(e);
    return t.copy(e), n > this.radius * this.radius && (t.sub(this.center).normalize(), t.multiplyScalar(this.radius).add(this.center)), t;
  }
  getBoundingBox(e) {
    return this.isEmpty() ? (e.makeEmpty(), e) : (e.set(this.center, this.center), e.expandByScalar(this.radius), e);
  }
  applyMatrix4(e) {
    return this.center.applyMatrix4(e), this.radius = this.radius * e.getMaxScaleOnAxis(), this;
  }
  translate(e) {
    return this.center.add(e), this;
  }
  expandByPoint(e) {
    if (this.isEmpty())
      return this.center.copy(e), this.radius = 0, this;
    oi.subVectors(e, this.center);
    const t = oi.lengthSq();
    if (t > this.radius * this.radius) {
      const n = Math.sqrt(t), r = (n - this.radius) * 0.5;
      this.center.addScaledVector(oi, r / n), this.radius += r;
    }
    return this;
  }
  union(e) {
    return e.isEmpty() ? this : this.isEmpty() ? (this.copy(e), this) : (this.center.equals(e.center) === !0 ? this.radius = Math.max(this.radius, e.radius) : (Ar.subVectors(e.center, this.center).setLength(e.radius), this.expandByPoint(oi.copy(e.center).add(Ar)), this.expandByPoint(oi.copy(e.center).sub(Ar))), this);
  }
  equals(e) {
    return e.center.equals(this.center) && e.radius === this.radius;
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
const Vt = /* @__PURE__ */ new U(), br = /* @__PURE__ */ new U(), Li = /* @__PURE__ */ new U(), $t = /* @__PURE__ */ new U(), wr = /* @__PURE__ */ new U(), Ui = /* @__PURE__ */ new U(), Rr = /* @__PURE__ */ new U();
class cs {
  constructor(e = new U(), t = new U(0, 0, -1)) {
    this.origin = e, this.direction = t;
  }
  set(e, t) {
    return this.origin.copy(e), this.direction.copy(t), this;
  }
  copy(e) {
    return this.origin.copy(e.origin), this.direction.copy(e.direction), this;
  }
  at(e, t) {
    return t.copy(this.origin).addScaledVector(this.direction, e);
  }
  lookAt(e) {
    return this.direction.copy(e).sub(this.origin).normalize(), this;
  }
  recast(e) {
    return this.origin.copy(this.at(e, Vt)), this;
  }
  closestPointToPoint(e, t) {
    t.subVectors(e, this.origin);
    const n = t.dot(this.direction);
    return n < 0 ? t.copy(this.origin) : t.copy(this.origin).addScaledVector(this.direction, n);
  }
  distanceToPoint(e) {
    return Math.sqrt(this.distanceSqToPoint(e));
  }
  distanceSqToPoint(e) {
    const t = Vt.subVectors(e, this.origin).dot(this.direction);
    return t < 0 ? this.origin.distanceToSquared(e) : (Vt.copy(this.origin).addScaledVector(this.direction, t), Vt.distanceToSquared(e));
  }
  distanceSqToSegment(e, t, n, r) {
    br.copy(e).add(t).multiplyScalar(0.5), Li.copy(t).sub(e).normalize(), $t.copy(this.origin).sub(br);
    const s = e.distanceTo(t) * 0.5, o = -this.direction.dot(Li), a = $t.dot(this.direction), l = -$t.dot(Li), c = $t.lengthSq(), h = Math.abs(1 - o * o);
    let f, u, m, g;
    if (h > 0)
      if (f = o * l - a, u = o * a - l, g = s * h, f >= 0)
        if (u >= -g)
          if (u <= g) {
            const x = 1 / h;
            f *= x, u *= x, m = f * (f + o * u + 2 * a) + u * (o * f + u + 2 * l) + c;
          } else
            u = s, f = Math.max(0, -(o * u + a)), m = -f * f + u * (u + 2 * l) + c;
        else
          u = -s, f = Math.max(0, -(o * u + a)), m = -f * f + u * (u + 2 * l) + c;
      else
        u <= -g ? (f = Math.max(0, -(-o * s + a)), u = f > 0 ? -s : Math.min(Math.max(-s, -l), s), m = -f * f + u * (u + 2 * l) + c) : u <= g ? (f = 0, u = Math.min(Math.max(-s, -l), s), m = u * (u + 2 * l) + c) : (f = Math.max(0, -(o * s + a)), u = f > 0 ? s : Math.min(Math.max(-s, -l), s), m = -f * f + u * (u + 2 * l) + c);
    else
      u = o > 0 ? -s : s, f = Math.max(0, -(o * u + a)), m = -f * f + u * (u + 2 * l) + c;
    return n && n.copy(this.origin).addScaledVector(this.direction, f), r && r.copy(br).addScaledVector(Li, u), m;
  }
  intersectSphere(e, t) {
    Vt.subVectors(e.center, this.origin);
    const n = Vt.dot(this.direction), r = Vt.dot(Vt) - n * n, s = e.radius * e.radius;
    if (r > s)
      return null;
    const o = Math.sqrt(s - r), a = n - o, l = n + o;
    return l < 0 ? null : a < 0 ? this.at(l, t) : this.at(a, t);
  }
  intersectsSphere(e) {
    return this.distanceSqToPoint(e.center) <= e.radius * e.radius;
  }
  distanceToPlane(e) {
    const t = e.normal.dot(this.direction);
    if (t === 0)
      return e.distanceToPoint(this.origin) === 0 ? 0 : null;
    const n = -(this.origin.dot(e.normal) + e.constant) / t;
    return n >= 0 ? n : null;
  }
  intersectPlane(e, t) {
    const n = this.distanceToPlane(e);
    return n === null ? null : this.at(n, t);
  }
  intersectsPlane(e) {
    const t = e.distanceToPoint(this.origin);
    return t === 0 || e.normal.dot(this.direction) * t < 0;
  }
  intersectBox(e, t) {
    let n, r, s, o, a, l;
    const c = 1 / this.direction.x, h = 1 / this.direction.y, f = 1 / this.direction.z, u = this.origin;
    return c >= 0 ? (n = (e.min.x - u.x) * c, r = (e.max.x - u.x) * c) : (n = (e.max.x - u.x) * c, r = (e.min.x - u.x) * c), h >= 0 ? (s = (e.min.y - u.y) * h, o = (e.max.y - u.y) * h) : (s = (e.max.y - u.y) * h, o = (e.min.y - u.y) * h), n > o || s > r || ((s > n || isNaN(n)) && (n = s), (o < r || isNaN(r)) && (r = o), f >= 0 ? (a = (e.min.z - u.z) * f, l = (e.max.z - u.z) * f) : (a = (e.max.z - u.z) * f, l = (e.min.z - u.z) * f), n > l || a > r) || ((a > n || n !== n) && (n = a), (l < r || r !== r) && (r = l), r < 0) ? null : this.at(n >= 0 ? n : r, t);
  }
  intersectsBox(e) {
    return this.intersectBox(e, Vt) !== null;
  }
  intersectTriangle(e, t, n, r, s) {
    wr.subVectors(t, e), Ui.subVectors(n, e), Rr.crossVectors(wr, Ui);
    let o = this.direction.dot(Rr), a;
    if (o > 0) {
      if (r)
        return null;
      a = 1;
    } else if (o < 0)
      a = -1, o = -o;
    else
      return null;
    $t.subVectors(this.origin, e);
    const l = a * this.direction.dot(Ui.crossVectors($t, Ui));
    if (l < 0)
      return null;
    const c = a * this.direction.dot(wr.cross($t));
    if (c < 0 || l + c > o)
      return null;
    const h = -a * $t.dot(Rr);
    return h < 0 ? null : this.at(h / o, s);
  }
  applyMatrix4(e) {
    return this.origin.applyMatrix4(e), this.direction.transformDirection(e), this;
  }
  equals(e) {
    return e.origin.equals(this.origin) && e.direction.equals(this.direction);
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
class nt {
  constructor(e, t, n, r, s, o, a, l, c, h, f, u, m, g, x, p) {
    nt.prototype.isMatrix4 = !0, this.elements = [
      1,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1
    ], e !== void 0 && this.set(e, t, n, r, s, o, a, l, c, h, f, u, m, g, x, p);
  }
  set(e, t, n, r, s, o, a, l, c, h, f, u, m, g, x, p) {
    const d = this.elements;
    return d[0] = e, d[4] = t, d[8] = n, d[12] = r, d[1] = s, d[5] = o, d[9] = a, d[13] = l, d[2] = c, d[6] = h, d[10] = f, d[14] = u, d[3] = m, d[7] = g, d[11] = x, d[15] = p, this;
  }
  identity() {
    return this.set(
      1,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  clone() {
    return new nt().fromArray(this.elements);
  }
  copy(e) {
    const t = this.elements, n = e.elements;
    return t[0] = n[0], t[1] = n[1], t[2] = n[2], t[3] = n[3], t[4] = n[4], t[5] = n[5], t[6] = n[6], t[7] = n[7], t[8] = n[8], t[9] = n[9], t[10] = n[10], t[11] = n[11], t[12] = n[12], t[13] = n[13], t[14] = n[14], t[15] = n[15], this;
  }
  copyPosition(e) {
    const t = this.elements, n = e.elements;
    return t[12] = n[12], t[13] = n[13], t[14] = n[14], this;
  }
  setFromMatrix3(e) {
    const t = e.elements;
    return this.set(
      t[0],
      t[3],
      t[6],
      0,
      t[1],
      t[4],
      t[7],
      0,
      t[2],
      t[5],
      t[8],
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  extractBasis(e, t, n) {
    return e.setFromMatrixColumn(this, 0), t.setFromMatrixColumn(this, 1), n.setFromMatrixColumn(this, 2), this;
  }
  makeBasis(e, t, n) {
    return this.set(
      e.x,
      t.x,
      n.x,
      0,
      e.y,
      t.y,
      n.y,
      0,
      e.z,
      t.z,
      n.z,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  extractRotation(e) {
    const t = this.elements, n = e.elements, r = 1 / Fn.setFromMatrixColumn(e, 0).length(), s = 1 / Fn.setFromMatrixColumn(e, 1).length(), o = 1 / Fn.setFromMatrixColumn(e, 2).length();
    return t[0] = n[0] * r, t[1] = n[1] * r, t[2] = n[2] * r, t[3] = 0, t[4] = n[4] * s, t[5] = n[5] * s, t[6] = n[6] * s, t[7] = 0, t[8] = n[8] * o, t[9] = n[9] * o, t[10] = n[10] * o, t[11] = 0, t[12] = 0, t[13] = 0, t[14] = 0, t[15] = 1, this;
  }
  makeRotationFromEuler(e) {
    const t = this.elements, n = e.x, r = e.y, s = e.z, o = Math.cos(n), a = Math.sin(n), l = Math.cos(r), c = Math.sin(r), h = Math.cos(s), f = Math.sin(s);
    if (e.order === "XYZ") {
      const u = o * h, m = o * f, g = a * h, x = a * f;
      t[0] = l * h, t[4] = -l * f, t[8] = c, t[1] = m + g * c, t[5] = u - x * c, t[9] = -a * l, t[2] = x - u * c, t[6] = g + m * c, t[10] = o * l;
    } else if (e.order === "YXZ") {
      const u = l * h, m = l * f, g = c * h, x = c * f;
      t[0] = u + x * a, t[4] = g * a - m, t[8] = o * c, t[1] = o * f, t[5] = o * h, t[9] = -a, t[2] = m * a - g, t[6] = x + u * a, t[10] = o * l;
    } else if (e.order === "ZXY") {
      const u = l * h, m = l * f, g = c * h, x = c * f;
      t[0] = u - x * a, t[4] = -o * f, t[8] = g + m * a, t[1] = m + g * a, t[5] = o * h, t[9] = x - u * a, t[2] = -o * c, t[6] = a, t[10] = o * l;
    } else if (e.order === "ZYX") {
      const u = o * h, m = o * f, g = a * h, x = a * f;
      t[0] = l * h, t[4] = g * c - m, t[8] = u * c + x, t[1] = l * f, t[5] = x * c + u, t[9] = m * c - g, t[2] = -c, t[6] = a * l, t[10] = o * l;
    } else if (e.order === "YZX") {
      const u = o * l, m = o * c, g = a * l, x = a * c;
      t[0] = l * h, t[4] = x - u * f, t[8] = g * f + m, t[1] = f, t[5] = o * h, t[9] = -a * h, t[2] = -c * h, t[6] = m * f + g, t[10] = u - x * f;
    } else if (e.order === "XZY") {
      const u = o * l, m = o * c, g = a * l, x = a * c;
      t[0] = l * h, t[4] = -f, t[8] = c * h, t[1] = u * f + x, t[5] = o * h, t[9] = m * f - g, t[2] = g * f - m, t[6] = a * h, t[10] = x * f + u;
    }
    return t[3] = 0, t[7] = 0, t[11] = 0, t[12] = 0, t[13] = 0, t[14] = 0, t[15] = 1, this;
  }
  makeRotationFromQuaternion(e) {
    return this.compose(vc, e, xc);
  }
  lookAt(e, t, n) {
    const r = this.elements;
    return vt.subVectors(e, t), vt.lengthSq() === 0 && (vt.z = 1), vt.normalize(), Qt.crossVectors(n, vt), Qt.lengthSq() === 0 && (Math.abs(n.z) === 1 ? vt.x += 1e-4 : vt.z += 1e-4, vt.normalize(), Qt.crossVectors(n, vt)), Qt.normalize(), Di.crossVectors(vt, Qt), r[0] = Qt.x, r[4] = Di.x, r[8] = vt.x, r[1] = Qt.y, r[5] = Di.y, r[9] = vt.y, r[2] = Qt.z, r[6] = Di.z, r[10] = vt.z, this;
  }
  multiply(e) {
    return this.multiplyMatrices(this, e);
  }
  premultiply(e) {
    return this.multiplyMatrices(e, this);
  }
  multiplyMatrices(e, t) {
    const n = e.elements, r = t.elements, s = this.elements, o = n[0], a = n[4], l = n[8], c = n[12], h = n[1], f = n[5], u = n[9], m = n[13], g = n[2], x = n[6], p = n[10], d = n[14], A = n[3], _ = n[7], T = n[11], C = n[15], L = r[0], w = r[4], V = r[8], M = r[12], y = r[1], Y = r[5], ce = r[9], B = r[13], H = r[2], G = r[6], Q = r[10], X = r[14], j = r[3], J = r[7], ee = r[11], I = r[15];
    return s[0] = o * L + a * y + l * H + c * j, s[4] = o * w + a * Y + l * G + c * J, s[8] = o * V + a * ce + l * Q + c * ee, s[12] = o * M + a * B + l * X + c * I, s[1] = h * L + f * y + u * H + m * j, s[5] = h * w + f * Y + u * G + m * J, s[9] = h * V + f * ce + u * Q + m * ee, s[13] = h * M + f * B + u * X + m * I, s[2] = g * L + x * y + p * H + d * j, s[6] = g * w + x * Y + p * G + d * J, s[10] = g * V + x * ce + p * Q + d * ee, s[14] = g * M + x * B + p * X + d * I, s[3] = A * L + _ * y + T * H + C * j, s[7] = A * w + _ * Y + T * G + C * J, s[11] = A * V + _ * ce + T * Q + C * ee, s[15] = A * M + _ * B + T * X + C * I, this;
  }
  multiplyScalar(e) {
    const t = this.elements;
    return t[0] *= e, t[4] *= e, t[8] *= e, t[12] *= e, t[1] *= e, t[5] *= e, t[9] *= e, t[13] *= e, t[2] *= e, t[6] *= e, t[10] *= e, t[14] *= e, t[3] *= e, t[7] *= e, t[11] *= e, t[15] *= e, this;
  }
  determinant() {
    const e = this.elements, t = e[0], n = e[4], r = e[8], s = e[12], o = e[1], a = e[5], l = e[9], c = e[13], h = e[2], f = e[6], u = e[10], m = e[14], g = e[3], x = e[7], p = e[11], d = e[15];
    return g * (+s * l * f - r * c * f - s * a * u + n * c * u + r * a * m - n * l * m) + x * (+t * l * m - t * c * u + s * o * u - r * o * m + r * c * h - s * l * h) + p * (+t * c * f - t * a * m - s * o * f + n * o * m + s * a * h - n * c * h) + d * (-r * a * h - t * l * f + t * a * u + r * o * f - n * o * u + n * l * h);
  }
  transpose() {
    const e = this.elements;
    let t;
    return t = e[1], e[1] = e[4], e[4] = t, t = e[2], e[2] = e[8], e[8] = t, t = e[6], e[6] = e[9], e[9] = t, t = e[3], e[3] = e[12], e[12] = t, t = e[7], e[7] = e[13], e[13] = t, t = e[11], e[11] = e[14], e[14] = t, this;
  }
  setPosition(e, t, n) {
    const r = this.elements;
    return e.isVector3 ? (r[12] = e.x, r[13] = e.y, r[14] = e.z) : (r[12] = e, r[13] = t, r[14] = n), this;
  }
  invert() {
    const e = this.elements, t = e[0], n = e[1], r = e[2], s = e[3], o = e[4], a = e[5], l = e[6], c = e[7], h = e[8], f = e[9], u = e[10], m = e[11], g = e[12], x = e[13], p = e[14], d = e[15], A = f * p * c - x * u * c + x * l * m - a * p * m - f * l * d + a * u * d, _ = g * u * c - h * p * c - g * l * m + o * p * m + h * l * d - o * u * d, T = h * x * c - g * f * c + g * a * m - o * x * m - h * a * d + o * f * d, C = g * f * l - h * x * l - g * a * u + o * x * u + h * a * p - o * f * p, L = t * A + n * _ + r * T + s * C;
    if (L === 0)
      return this.set(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    const w = 1 / L;
    return e[0] = A * w, e[1] = (x * u * s - f * p * s - x * r * m + n * p * m + f * r * d - n * u * d) * w, e[2] = (a * p * s - x * l * s + x * r * c - n * p * c - a * r * d + n * l * d) * w, e[3] = (f * l * s - a * u * s - f * r * c + n * u * c + a * r * m - n * l * m) * w, e[4] = _ * w, e[5] = (h * p * s - g * u * s + g * r * m - t * p * m - h * r * d + t * u * d) * w, e[6] = (g * l * s - o * p * s - g * r * c + t * p * c + o * r * d - t * l * d) * w, e[7] = (o * u * s - h * l * s + h * r * c - t * u * c - o * r * m + t * l * m) * w, e[8] = T * w, e[9] = (g * f * s - h * x * s - g * n * m + t * x * m + h * n * d - t * f * d) * w, e[10] = (o * x * s - g * a * s + g * n * c - t * x * c - o * n * d + t * a * d) * w, e[11] = (h * a * s - o * f * s - h * n * c + t * f * c + o * n * m - t * a * m) * w, e[12] = C * w, e[13] = (h * x * r - g * f * r + g * n * u - t * x * u - h * n * p + t * f * p) * w, e[14] = (g * a * r - o * x * r - g * n * l + t * x * l + o * n * p - t * a * p) * w, e[15] = (o * f * r - h * a * r + h * n * l - t * f * l - o * n * u + t * a * u) * w, this;
  }
  scale(e) {
    const t = this.elements, n = e.x, r = e.y, s = e.z;
    return t[0] *= n, t[4] *= r, t[8] *= s, t[1] *= n, t[5] *= r, t[9] *= s, t[2] *= n, t[6] *= r, t[10] *= s, t[3] *= n, t[7] *= r, t[11] *= s, this;
  }
  getMaxScaleOnAxis() {
    const e = this.elements, t = e[0] * e[0] + e[1] * e[1] + e[2] * e[2], n = e[4] * e[4] + e[5] * e[5] + e[6] * e[6], r = e[8] * e[8] + e[9] * e[9] + e[10] * e[10];
    return Math.sqrt(Math.max(t, n, r));
  }
  makeTranslation(e, t, n) {
    return e.isVector3 ? this.set(
      1,
      0,
      0,
      e.x,
      0,
      1,
      0,
      e.y,
      0,
      0,
      1,
      e.z,
      0,
      0,
      0,
      1
    ) : this.set(
      1,
      0,
      0,
      e,
      0,
      1,
      0,
      t,
      0,
      0,
      1,
      n,
      0,
      0,
      0,
      1
    ), this;
  }
  makeRotationX(e) {
    const t = Math.cos(e), n = Math.sin(e);
    return this.set(
      1,
      0,
      0,
      0,
      0,
      t,
      -n,
      0,
      0,
      n,
      t,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  makeRotationY(e) {
    const t = Math.cos(e), n = Math.sin(e);
    return this.set(
      t,
      0,
      n,
      0,
      0,
      1,
      0,
      0,
      -n,
      0,
      t,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  makeRotationZ(e) {
    const t = Math.cos(e), n = Math.sin(e);
    return this.set(
      t,
      -n,
      0,
      0,
      n,
      t,
      0,
      0,
      0,
      0,
      1,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  makeRotationAxis(e, t) {
    const n = Math.cos(t), r = Math.sin(t), s = 1 - n, o = e.x, a = e.y, l = e.z, c = s * o, h = s * a;
    return this.set(
      c * o + n,
      c * a - r * l,
      c * l + r * a,
      0,
      c * a + r * l,
      h * a + n,
      h * l - r * o,
      0,
      c * l - r * a,
      h * l + r * o,
      s * l * l + n,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  makeScale(e, t, n) {
    return this.set(
      e,
      0,
      0,
      0,
      0,
      t,
      0,
      0,
      0,
      0,
      n,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  makeShear(e, t, n, r, s, o) {
    return this.set(
      1,
      n,
      s,
      0,
      e,
      1,
      o,
      0,
      t,
      r,
      1,
      0,
      0,
      0,
      0,
      1
    ), this;
  }
  compose(e, t, n) {
    const r = this.elements, s = t._x, o = t._y, a = t._z, l = t._w, c = s + s, h = o + o, f = a + a, u = s * c, m = s * h, g = s * f, x = o * h, p = o * f, d = a * f, A = l * c, _ = l * h, T = l * f, C = n.x, L = n.y, w = n.z;
    return r[0] = (1 - (x + d)) * C, r[1] = (m + T) * C, r[2] = (g - _) * C, r[3] = 0, r[4] = (m - T) * L, r[5] = (1 - (u + d)) * L, r[6] = (p + A) * L, r[7] = 0, r[8] = (g + _) * w, r[9] = (p - A) * w, r[10] = (1 - (u + x)) * w, r[11] = 0, r[12] = e.x, r[13] = e.y, r[14] = e.z, r[15] = 1, this;
  }
  decompose(e, t, n) {
    const r = this.elements;
    let s = Fn.set(r[0], r[1], r[2]).length();
    const o = Fn.set(r[4], r[5], r[6]).length(), a = Fn.set(r[8], r[9], r[10]).length();
    this.determinant() < 0 && (s = -s), e.x = r[12], e.y = r[13], e.z = r[14], wt.copy(this);
    const c = 1 / s, h = 1 / o, f = 1 / a;
    return wt.elements[0] *= c, wt.elements[1] *= c, wt.elements[2] *= c, wt.elements[4] *= h, wt.elements[5] *= h, wt.elements[6] *= h, wt.elements[8] *= f, wt.elements[9] *= f, wt.elements[10] *= f, t.setFromRotationMatrix(wt), n.x = s, n.y = o, n.z = a, this;
  }
  makePerspective(e, t, n, r, s, o, a = qt) {
    const l = this.elements, c = 2 * s / (t - e), h = 2 * s / (n - r), f = (t + e) / (t - e), u = (n + r) / (n - r);
    let m, g;
    if (a === qt)
      m = -(o + s) / (o - s), g = -2 * o * s / (o - s);
    else if (a === tr)
      m = -o / (o - s), g = -o * s / (o - s);
    else
      throw new Error("THREE.Matrix4.makePerspective(): Invalid coordinate system: " + a);
    return l[0] = c, l[4] = 0, l[8] = f, l[12] = 0, l[1] = 0, l[5] = h, l[9] = u, l[13] = 0, l[2] = 0, l[6] = 0, l[10] = m, l[14] = g, l[3] = 0, l[7] = 0, l[11] = -1, l[15] = 0, this;
  }
  makeOrthographic(e, t, n, r, s, o, a = qt) {
    const l = this.elements, c = 1 / (t - e), h = 1 / (n - r), f = 1 / (o - s), u = (t + e) * c, m = (n + r) * h;
    let g, x;
    if (a === qt)
      g = (o + s) * f, x = -2 * f;
    else if (a === tr)
      g = s * f, x = -1 * f;
    else
      throw new Error("THREE.Matrix4.makeOrthographic(): Invalid coordinate system: " + a);
    return l[0] = 2 * c, l[4] = 0, l[8] = 0, l[12] = -u, l[1] = 0, l[5] = 2 * h, l[9] = 0, l[13] = -m, l[2] = 0, l[6] = 0, l[10] = x, l[14] = -g, l[3] = 0, l[7] = 0, l[11] = 0, l[15] = 1, this;
  }
  equals(e) {
    const t = this.elements, n = e.elements;
    for (let r = 0; r < 16; r++)
      if (t[r] !== n[r])
        return !1;
    return !0;
  }
  fromArray(e, t = 0) {
    for (let n = 0; n < 16; n++)
      this.elements[n] = e[n + t];
    return this;
  }
  toArray(e = [], t = 0) {
    const n = this.elements;
    return e[t] = n[0], e[t + 1] = n[1], e[t + 2] = n[2], e[t + 3] = n[3], e[t + 4] = n[4], e[t + 5] = n[5], e[t + 6] = n[6], e[t + 7] = n[7], e[t + 8] = n[8], e[t + 9] = n[9], e[t + 10] = n[10], e[t + 11] = n[11], e[t + 12] = n[12], e[t + 13] = n[13], e[t + 14] = n[14], e[t + 15] = n[15], e;
  }
}
const Fn = /* @__PURE__ */ new U(), wt = /* @__PURE__ */ new nt(), vc = /* @__PURE__ */ new U(0, 0, 0), xc = /* @__PURE__ */ new U(1, 1, 1), Qt = /* @__PURE__ */ new U(), Di = /* @__PURE__ */ new U(), vt = /* @__PURE__ */ new U(), ra = /* @__PURE__ */ new nt(), sa = /* @__PURE__ */ new An();
class sr {
  constructor(e = 0, t = 0, n = 0, r = sr.DEFAULT_ORDER) {
    this.isEuler = !0, this._x = e, this._y = t, this._z = n, this._order = r;
  }
  get x() {
    return this._x;
  }
  set x(e) {
    this._x = e, this._onChangeCallback();
  }
  get y() {
    return this._y;
  }
  set y(e) {
    this._y = e, this._onChangeCallback();
  }
  get z() {
    return this._z;
  }
  set z(e) {
    this._z = e, this._onChangeCallback();
  }
  get order() {
    return this._order;
  }
  set order(e) {
    this._order = e, this._onChangeCallback();
  }
  set(e, t, n, r = this._order) {
    return this._x = e, this._y = t, this._z = n, this._order = r, this._onChangeCallback(), this;
  }
  clone() {
    return new this.constructor(this._x, this._y, this._z, this._order);
  }
  copy(e) {
    return this._x = e._x, this._y = e._y, this._z = e._z, this._order = e._order, this._onChangeCallback(), this;
  }
  setFromRotationMatrix(e, t = this._order, n = !0) {
    const r = e.elements, s = r[0], o = r[4], a = r[8], l = r[1], c = r[5], h = r[9], f = r[2], u = r[6], m = r[10];
    switch (t) {
      case "XYZ":
        this._y = Math.asin(it(a, -1, 1)), Math.abs(a) < 0.9999999 ? (this._x = Math.atan2(-h, m), this._z = Math.atan2(-o, s)) : (this._x = Math.atan2(u, c), this._z = 0);
        break;
      case "YXZ":
        this._x = Math.asin(-it(h, -1, 1)), Math.abs(h) < 0.9999999 ? (this._y = Math.atan2(a, m), this._z = Math.atan2(l, c)) : (this._y = Math.atan2(-f, s), this._z = 0);
        break;
      case "ZXY":
        this._x = Math.asin(it(u, -1, 1)), Math.abs(u) < 0.9999999 ? (this._y = Math.atan2(-f, m), this._z = Math.atan2(-o, c)) : (this._y = 0, this._z = Math.atan2(l, s));
        break;
      case "ZYX":
        this._y = Math.asin(-it(f, -1, 1)), Math.abs(f) < 0.9999999 ? (this._x = Math.atan2(u, m), this._z = Math.atan2(l, s)) : (this._x = 0, this._z = Math.atan2(-o, c));
        break;
      case "YZX":
        this._z = Math.asin(it(l, -1, 1)), Math.abs(l) < 0.9999999 ? (this._x = Math.atan2(-h, c), this._y = Math.atan2(-f, s)) : (this._x = 0, this._y = Math.atan2(a, m));
        break;
      case "XZY":
        this._z = Math.asin(-it(o, -1, 1)), Math.abs(o) < 0.9999999 ? (this._x = Math.atan2(u, c), this._y = Math.atan2(a, s)) : (this._x = Math.atan2(-h, m), this._y = 0);
        break;
      default:
        console.warn("THREE.Euler: .setFromRotationMatrix() encountered an unknown order: " + t);
    }
    return this._order = t, n === !0 && this._onChangeCallback(), this;
  }
  setFromQuaternion(e, t, n) {
    return ra.makeRotationFromQuaternion(e), this.setFromRotationMatrix(ra, t, n);
  }
  setFromVector3(e, t = this._order) {
    return this.set(e.x, e.y, e.z, t);
  }
  reorder(e) {
    return sa.setFromEuler(this), this.setFromQuaternion(sa, e);
  }
  equals(e) {
    return e._x === this._x && e._y === this._y && e._z === this._z && e._order === this._order;
  }
  fromArray(e) {
    return this._x = e[0], this._y = e[1], this._z = e[2], e[3] !== void 0 && (this._order = e[3]), this._onChangeCallback(), this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this._x, e[t + 1] = this._y, e[t + 2] = this._z, e[t + 3] = this._order, e;
  }
  _onChange(e) {
    return this._onChangeCallback = e, this;
  }
  _onChangeCallback() {
  }
  *[Symbol.iterator]() {
    yield this._x, yield this._y, yield this._z, yield this._order;
  }
}
sr.DEFAULT_ORDER = "XYZ";
class hs {
  constructor() {
    this.mask = 1;
  }
  set(e) {
    this.mask = (1 << e | 0) >>> 0;
  }
  enable(e) {
    this.mask |= 1 << e | 0;
  }
  enableAll() {
    this.mask = -1;
  }
  toggle(e) {
    this.mask ^= 1 << e | 0;
  }
  disable(e) {
    this.mask &= ~(1 << e | 0);
  }
  disableAll() {
    this.mask = 0;
  }
  test(e) {
    return (this.mask & e.mask) !== 0;
  }
  isEnabled(e) {
    return (this.mask & (1 << e | 0)) !== 0;
  }
}
let Mc = 0;
const aa = /* @__PURE__ */ new U(), Bn = /* @__PURE__ */ new An(), kt = /* @__PURE__ */ new nt(), Ii = /* @__PURE__ */ new U(), li = /* @__PURE__ */ new U(), Sc = /* @__PURE__ */ new U(), Ec = /* @__PURE__ */ new An(), oa = /* @__PURE__ */ new U(1, 0, 0), la = /* @__PURE__ */ new U(0, 1, 0), ca = /* @__PURE__ */ new U(0, 0, 1), yc = { type: "added" }, ha = { type: "removed" };
class ht extends Rn {
  constructor() {
    super(), this.isObject3D = !0, Object.defineProperty(this, "id", { value: Mc++ }), this.uuid = Cn(), this.name = "", this.type = "Object3D", this.parent = null, this.children = [], this.up = ht.DEFAULT_UP.clone();
    const e = new U(), t = new sr(), n = new An(), r = new U(1, 1, 1);
    function s() {
      n.setFromEuler(t, !1);
    }
    function o() {
      t.setFromQuaternion(n, void 0, !1);
    }
    t._onChange(s), n._onChange(o), Object.defineProperties(this, {
      position: {
        configurable: !0,
        enumerable: !0,
        value: e
      },
      rotation: {
        configurable: !0,
        enumerable: !0,
        value: t
      },
      quaternion: {
        configurable: !0,
        enumerable: !0,
        value: n
      },
      scale: {
        configurable: !0,
        enumerable: !0,
        value: r
      },
      modelViewMatrix: {
        value: new nt()
      },
      normalMatrix: {
        value: new Be()
      }
    }), this.matrix = new nt(), this.matrixWorld = new nt(), this.matrixAutoUpdate = ht.DEFAULT_MATRIX_AUTO_UPDATE, this.matrixWorldNeedsUpdate = !1, this.matrixWorldAutoUpdate = ht.DEFAULT_MATRIX_WORLD_AUTO_UPDATE, this.layers = new hs(), this.visible = !0, this.castShadow = !1, this.receiveShadow = !1, this.frustumCulled = !0, this.renderOrder = 0, this.animations = [], this.userData = {};
  }
  onBeforeRender() {
  }
  onAfterRender() {
  }
  applyMatrix4(e) {
    this.matrixAutoUpdate && this.updateMatrix(), this.matrix.premultiply(e), this.matrix.decompose(this.position, this.quaternion, this.scale);
  }
  applyQuaternion(e) {
    return this.quaternion.premultiply(e), this;
  }
  setRotationFromAxisAngle(e, t) {
    this.quaternion.setFromAxisAngle(e, t);
  }
  setRotationFromEuler(e) {
    this.quaternion.setFromEuler(e, !0);
  }
  setRotationFromMatrix(e) {
    this.quaternion.setFromRotationMatrix(e);
  }
  setRotationFromQuaternion(e) {
    this.quaternion.copy(e);
  }
  rotateOnAxis(e, t) {
    return Bn.setFromAxisAngle(e, t), this.quaternion.multiply(Bn), this;
  }
  rotateOnWorldAxis(e, t) {
    return Bn.setFromAxisAngle(e, t), this.quaternion.premultiply(Bn), this;
  }
  rotateX(e) {
    return this.rotateOnAxis(oa, e);
  }
  rotateY(e) {
    return this.rotateOnAxis(la, e);
  }
  rotateZ(e) {
    return this.rotateOnAxis(ca, e);
  }
  translateOnAxis(e, t) {
    return aa.copy(e).applyQuaternion(this.quaternion), this.position.add(aa.multiplyScalar(t)), this;
  }
  translateX(e) {
    return this.translateOnAxis(oa, e);
  }
  translateY(e) {
    return this.translateOnAxis(la, e);
  }
  translateZ(e) {
    return this.translateOnAxis(ca, e);
  }
  localToWorld(e) {
    return this.updateWorldMatrix(!0, !1), e.applyMatrix4(this.matrixWorld);
  }
  worldToLocal(e) {
    return this.updateWorldMatrix(!0, !1), e.applyMatrix4(kt.copy(this.matrixWorld).invert());
  }
  lookAt(e, t, n) {
    e.isVector3 ? Ii.copy(e) : Ii.set(e, t, n);
    const r = this.parent;
    this.updateWorldMatrix(!0, !1), li.setFromMatrixPosition(this.matrixWorld), this.isCamera || this.isLight ? kt.lookAt(li, Ii, this.up) : kt.lookAt(Ii, li, this.up), this.quaternion.setFromRotationMatrix(kt), r && (kt.extractRotation(r.matrixWorld), Bn.setFromRotationMatrix(kt), this.quaternion.premultiply(Bn.invert()));
  }
  add(e) {
    if (arguments.length > 1) {
      for (let t = 0; t < arguments.length; t++)
        this.add(arguments[t]);
      return this;
    }
    return e === this ? (console.error("THREE.Object3D.add: object can't be added as a child of itself.", e), this) : (e && e.isObject3D ? (e.parent !== null && e.parent.remove(e), e.parent = this, this.children.push(e), e.dispatchEvent(yc)) : console.error("THREE.Object3D.add: object not an instance of THREE.Object3D.", e), this);
  }
  remove(e) {
    if (arguments.length > 1) {
      for (let n = 0; n < arguments.length; n++)
        this.remove(arguments[n]);
      return this;
    }
    const t = this.children.indexOf(e);
    return t !== -1 && (e.parent = null, this.children.splice(t, 1), e.dispatchEvent(ha)), this;
  }
  removeFromParent() {
    const e = this.parent;
    return e !== null && e.remove(this), this;
  }
  clear() {
    for (let e = 0; e < this.children.length; e++) {
      const t = this.children[e];
      t.parent = null, t.dispatchEvent(ha);
    }
    return this.children.length = 0, this;
  }
  attach(e) {
    return this.updateWorldMatrix(!0, !1), kt.copy(this.matrixWorld).invert(), e.parent !== null && (e.parent.updateWorldMatrix(!0, !1), kt.multiply(e.parent.matrixWorld)), e.applyMatrix4(kt), this.add(e), e.updateWorldMatrix(!1, !0), this;
  }
  getObjectById(e) {
    return this.getObjectByProperty("id", e);
  }
  getObjectByName(e) {
    return this.getObjectByProperty("name", e);
  }
  getObjectByProperty(e, t) {
    if (this[e] === t)
      return this;
    for (let n = 0, r = this.children.length; n < r; n++) {
      const o = this.children[n].getObjectByProperty(e, t);
      if (o !== void 0)
        return o;
    }
  }
  getObjectsByProperty(e, t) {
    let n = [];
    this[e] === t && n.push(this);
    for (let r = 0, s = this.children.length; r < s; r++) {
      const o = this.children[r].getObjectsByProperty(e, t);
      o.length > 0 && (n = n.concat(o));
    }
    return n;
  }
  getWorldPosition(e) {
    return this.updateWorldMatrix(!0, !1), e.setFromMatrixPosition(this.matrixWorld);
  }
  getWorldQuaternion(e) {
    return this.updateWorldMatrix(!0, !1), this.matrixWorld.decompose(li, e, Sc), e;
  }
  getWorldScale(e) {
    return this.updateWorldMatrix(!0, !1), this.matrixWorld.decompose(li, Ec, e), e;
  }
  getWorldDirection(e) {
    this.updateWorldMatrix(!0, !1);
    const t = this.matrixWorld.elements;
    return e.set(t[8], t[9], t[10]).normalize();
  }
  raycast() {
  }
  traverse(e) {
    e(this);
    const t = this.children;
    for (let n = 0, r = t.length; n < r; n++)
      t[n].traverse(e);
  }
  traverseVisible(e) {
    if (this.visible === !1)
      return;
    e(this);
    const t = this.children;
    for (let n = 0, r = t.length; n < r; n++)
      t[n].traverseVisible(e);
  }
  traverseAncestors(e) {
    const t = this.parent;
    t !== null && (e(t), t.traverseAncestors(e));
  }
  updateMatrix() {
    this.matrix.compose(this.position, this.quaternion, this.scale), this.matrixWorldNeedsUpdate = !0;
  }
  updateMatrixWorld(e) {
    this.matrixAutoUpdate && this.updateMatrix(), (this.matrixWorldNeedsUpdate || e) && (this.parent === null ? this.matrixWorld.copy(this.matrix) : this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix), this.matrixWorldNeedsUpdate = !1, e = !0);
    const t = this.children;
    for (let n = 0, r = t.length; n < r; n++) {
      const s = t[n];
      (s.matrixWorldAutoUpdate === !0 || e === !0) && s.updateMatrixWorld(e);
    }
  }
  updateWorldMatrix(e, t) {
    const n = this.parent;
    if (e === !0 && n !== null && n.matrixWorldAutoUpdate === !0 && n.updateWorldMatrix(!0, !1), this.matrixAutoUpdate && this.updateMatrix(), this.parent === null ? this.matrixWorld.copy(this.matrix) : this.matrixWorld.multiplyMatrices(this.parent.matrixWorld, this.matrix), t === !0) {
      const r = this.children;
      for (let s = 0, o = r.length; s < o; s++) {
        const a = r[s];
        a.matrixWorldAutoUpdate === !0 && a.updateWorldMatrix(!1, !0);
      }
    }
  }
  toJSON(e) {
    const t = e === void 0 || typeof e == "string", n = {};
    t && (e = {
      geometries: {},
      materials: {},
      textures: {},
      images: {},
      shapes: {},
      skeletons: {},
      animations: {},
      nodes: {}
    }, n.metadata = {
      version: 4.6,
      type: "Object",
      generator: "Object3D.toJSON"
    });
    const r = {};
    r.uuid = this.uuid, r.type = this.type, this.name !== "" && (r.name = this.name), this.castShadow === !0 && (r.castShadow = !0), this.receiveShadow === !0 && (r.receiveShadow = !0), this.visible === !1 && (r.visible = !1), this.frustumCulled === !1 && (r.frustumCulled = !1), this.renderOrder !== 0 && (r.renderOrder = this.renderOrder), Object.keys(this.userData).length > 0 && (r.userData = this.userData), r.layers = this.layers.mask, r.matrix = this.matrix.toArray(), r.up = this.up.toArray(), this.matrixAutoUpdate === !1 && (r.matrixAutoUpdate = !1), this.isInstancedMesh && (r.type = "InstancedMesh", r.count = this.count, r.instanceMatrix = this.instanceMatrix.toJSON(), this.instanceColor !== null && (r.instanceColor = this.instanceColor.toJSON()));
    function s(a, l) {
      return a[l.uuid] === void 0 && (a[l.uuid] = l.toJSON(e)), l.uuid;
    }
    if (this.isScene)
      this.background && (this.background.isColor ? r.background = this.background.toJSON() : this.background.isTexture && (r.background = this.background.toJSON(e).uuid)), this.environment && this.environment.isTexture && this.environment.isRenderTargetTexture !== !0 && (r.environment = this.environment.toJSON(e).uuid);
    else if (this.isMesh || this.isLine || this.isPoints) {
      r.geometry = s(e.geometries, this.geometry);
      const a = this.geometry.parameters;
      if (a !== void 0 && a.shapes !== void 0) {
        const l = a.shapes;
        if (Array.isArray(l))
          for (let c = 0, h = l.length; c < h; c++) {
            const f = l[c];
            s(e.shapes, f);
          }
        else
          s(e.shapes, l);
      }
    }
    if (this.isSkinnedMesh && (r.bindMode = this.bindMode, r.bindMatrix = this.bindMatrix.toArray(), this.skeleton !== void 0 && (s(e.skeletons, this.skeleton), r.skeleton = this.skeleton.uuid)), this.material !== void 0)
      if (Array.isArray(this.material)) {
        const a = [];
        for (let l = 0, c = this.material.length; l < c; l++)
          a.push(s(e.materials, this.material[l]));
        r.material = a;
      } else
        r.material = s(e.materials, this.material);
    if (this.children.length > 0) {
      r.children = [];
      for (let a = 0; a < this.children.length; a++)
        r.children.push(this.children[a].toJSON(e).object);
    }
    if (this.animations.length > 0) {
      r.animations = [];
      for (let a = 0; a < this.animations.length; a++) {
        const l = this.animations[a];
        r.animations.push(s(e.animations, l));
      }
    }
    if (t) {
      const a = o(e.geometries), l = o(e.materials), c = o(e.textures), h = o(e.images), f = o(e.shapes), u = o(e.skeletons), m = o(e.animations), g = o(e.nodes);
      a.length > 0 && (n.geometries = a), l.length > 0 && (n.materials = l), c.length > 0 && (n.textures = c), h.length > 0 && (n.images = h), f.length > 0 && (n.shapes = f), u.length > 0 && (n.skeletons = u), m.length > 0 && (n.animations = m), g.length > 0 && (n.nodes = g);
    }
    return n.object = r, n;
    function o(a) {
      const l = [];
      for (const c in a) {
        const h = a[c];
        delete h.metadata, l.push(h);
      }
      return l;
    }
  }
  clone(e) {
    return new this.constructor().copy(this, e);
  }
  copy(e, t = !0) {
    if (this.name = e.name, this.up.copy(e.up), this.position.copy(e.position), this.rotation.order = e.rotation.order, this.quaternion.copy(e.quaternion), this.scale.copy(e.scale), this.matrix.copy(e.matrix), this.matrixWorld.copy(e.matrixWorld), this.matrixAutoUpdate = e.matrixAutoUpdate, this.matrixWorldNeedsUpdate = e.matrixWorldNeedsUpdate, this.matrixWorldAutoUpdate = e.matrixWorldAutoUpdate, this.layers.mask = e.layers.mask, this.visible = e.visible, this.castShadow = e.castShadow, this.receiveShadow = e.receiveShadow, this.frustumCulled = e.frustumCulled, this.renderOrder = e.renderOrder, this.animations = e.animations.slice(), this.userData = JSON.parse(JSON.stringify(e.userData)), t === !0)
      for (let n = 0; n < e.children.length; n++) {
        const r = e.children[n];
        this.add(r.clone());
      }
    return this;
  }
}
ht.DEFAULT_UP = /* @__PURE__ */ new U(0, 1, 0);
ht.DEFAULT_MATRIX_AUTO_UPDATE = !0;
ht.DEFAULT_MATRIX_WORLD_AUTO_UPDATE = !0;
const Rt = /* @__PURE__ */ new U(), Wt = /* @__PURE__ */ new U(), Cr = /* @__PURE__ */ new U(), Xt = /* @__PURE__ */ new U(), zn = /* @__PURE__ */ new U(), Hn = /* @__PURE__ */ new U(), ua = /* @__PURE__ */ new U(), Pr = /* @__PURE__ */ new U(), Lr = /* @__PURE__ */ new U(), Ur = /* @__PURE__ */ new U();
let Ni = !1;
class Pt {
  constructor(e = new U(), t = new U(), n = new U()) {
    this.a = e, this.b = t, this.c = n;
  }
  static getNormal(e, t, n, r) {
    r.subVectors(n, t), Rt.subVectors(e, t), r.cross(Rt);
    const s = r.lengthSq();
    return s > 0 ? r.multiplyScalar(1 / Math.sqrt(s)) : r.set(0, 0, 0);
  }
  // static/instance method to calculate barycentric coordinates
  // based on: http://www.blackpawn.com/texts/pointinpoly/default.html
  static getBarycoord(e, t, n, r, s) {
    Rt.subVectors(r, t), Wt.subVectors(n, t), Cr.subVectors(e, t);
    const o = Rt.dot(Rt), a = Rt.dot(Wt), l = Rt.dot(Cr), c = Wt.dot(Wt), h = Wt.dot(Cr), f = o * c - a * a;
    if (f === 0)
      return s.set(-2, -1, -1);
    const u = 1 / f, m = (c * l - a * h) * u, g = (o * h - a * l) * u;
    return s.set(1 - m - g, g, m);
  }
  static containsPoint(e, t, n, r) {
    return this.getBarycoord(e, t, n, r, Xt), Xt.x >= 0 && Xt.y >= 0 && Xt.x + Xt.y <= 1;
  }
  static getUV(e, t, n, r, s, o, a, l) {
    return Ni === !1 && (console.warn("THREE.Triangle.getUV() has been renamed to THREE.Triangle.getInterpolation()."), Ni = !0), this.getInterpolation(e, t, n, r, s, o, a, l);
  }
  static getInterpolation(e, t, n, r, s, o, a, l) {
    return this.getBarycoord(e, t, n, r, Xt), l.setScalar(0), l.addScaledVector(s, Xt.x), l.addScaledVector(o, Xt.y), l.addScaledVector(a, Xt.z), l;
  }
  static isFrontFacing(e, t, n, r) {
    return Rt.subVectors(n, t), Wt.subVectors(e, t), Rt.cross(Wt).dot(r) < 0;
  }
  set(e, t, n) {
    return this.a.copy(e), this.b.copy(t), this.c.copy(n), this;
  }
  setFromPointsAndIndices(e, t, n, r) {
    return this.a.copy(e[t]), this.b.copy(e[n]), this.c.copy(e[r]), this;
  }
  setFromAttributeAndIndices(e, t, n, r) {
    return this.a.fromBufferAttribute(e, t), this.b.fromBufferAttribute(e, n), this.c.fromBufferAttribute(e, r), this;
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    return this.a.copy(e.a), this.b.copy(e.b), this.c.copy(e.c), this;
  }
  getArea() {
    return Rt.subVectors(this.c, this.b), Wt.subVectors(this.a, this.b), Rt.cross(Wt).length() * 0.5;
  }
  getMidpoint(e) {
    return e.addVectors(this.a, this.b).add(this.c).multiplyScalar(1 / 3);
  }
  getNormal(e) {
    return Pt.getNormal(this.a, this.b, this.c, e);
  }
  getPlane(e) {
    return e.setFromCoplanarPoints(this.a, this.b, this.c);
  }
  getBarycoord(e, t) {
    return Pt.getBarycoord(e, this.a, this.b, this.c, t);
  }
  getUV(e, t, n, r, s) {
    return Ni === !1 && (console.warn("THREE.Triangle.getUV() has been renamed to THREE.Triangle.getInterpolation()."), Ni = !0), Pt.getInterpolation(e, this.a, this.b, this.c, t, n, r, s);
  }
  getInterpolation(e, t, n, r, s) {
    return Pt.getInterpolation(e, this.a, this.b, this.c, t, n, r, s);
  }
  containsPoint(e) {
    return Pt.containsPoint(e, this.a, this.b, this.c);
  }
  isFrontFacing(e) {
    return Pt.isFrontFacing(this.a, this.b, this.c, e);
  }
  intersectsBox(e) {
    return e.intersectsTriangle(this);
  }
  closestPointToPoint(e, t) {
    const n = this.a, r = this.b, s = this.c;
    let o, a;
    zn.subVectors(r, n), Hn.subVectors(s, n), Pr.subVectors(e, n);
    const l = zn.dot(Pr), c = Hn.dot(Pr);
    if (l <= 0 && c <= 0)
      return t.copy(n);
    Lr.subVectors(e, r);
    const h = zn.dot(Lr), f = Hn.dot(Lr);
    if (h >= 0 && f <= h)
      return t.copy(r);
    const u = l * f - h * c;
    if (u <= 0 && l >= 0 && h <= 0)
      return o = l / (l - h), t.copy(n).addScaledVector(zn, o);
    Ur.subVectors(e, s);
    const m = zn.dot(Ur), g = Hn.dot(Ur);
    if (g >= 0 && m <= g)
      return t.copy(s);
    const x = m * c - l * g;
    if (x <= 0 && c >= 0 && g <= 0)
      return a = c / (c - g), t.copy(n).addScaledVector(Hn, a);
    const p = h * g - m * f;
    if (p <= 0 && f - h >= 0 && m - g >= 0)
      return ua.subVectors(s, r), a = (f - h) / (f - h + (m - g)), t.copy(r).addScaledVector(ua, a);
    const d = 1 / (p + x + u);
    return o = x * d, a = u * d, t.copy(n).addScaledVector(zn, o).addScaledVector(Hn, a);
  }
  equals(e) {
    return e.a.equals(this.a) && e.b.equals(this.b) && e.c.equals(this.c);
  }
}
let Tc = 0;
class Ai extends Rn {
  constructor() {
    super(), this.isMaterial = !0, Object.defineProperty(this, "id", { value: Tc++ }), this.uuid = Cn(), this.name = "", this.type = "Material", this.blending = $n, this.side = ln, this.vertexColors = !1, this.opacity = 1, this.transparent = !1, this.alphaHash = !1, this.blendSrc = Qa, this.blendDst = eo, this.blendEquation = jn, this.blendSrcAlpha = null, this.blendDstAlpha = null, this.blendEquationAlpha = null, this.depthFunc = Yr, this.depthTest = !0, this.depthWrite = !0, this.stencilWriteMask = 255, this.stencilFunc = Fl, this.stencilRef = 0, this.stencilFuncMask = 255, this.stencilFail = xr, this.stencilZFail = xr, this.stencilZPass = xr, this.stencilWrite = !1, this.clippingPlanes = null, this.clipIntersection = !1, this.clipShadows = !1, this.shadowSide = null, this.colorWrite = !0, this.precision = null, this.polygonOffset = !1, this.polygonOffsetFactor = 0, this.polygonOffsetUnits = 0, this.dithering = !1, this.alphaToCoverage = !1, this.premultipliedAlpha = !1, this.forceSinglePass = !1, this.visible = !0, this.toneMapped = !0, this.userData = {}, this.version = 0, this._alphaTest = 0;
  }
  get alphaTest() {
    return this._alphaTest;
  }
  set alphaTest(e) {
    this._alphaTest > 0 != e > 0 && this.version++, this._alphaTest = e;
  }
  onBuild() {
  }
  onBeforeRender() {
  }
  onBeforeCompile() {
  }
  customProgramCacheKey() {
    return this.onBeforeCompile.toString();
  }
  setValues(e) {
    if (e !== void 0)
      for (const t in e) {
        const n = e[t];
        if (n === void 0) {
          console.warn(`THREE.Material: parameter '${t}' has value of undefined.`);
          continue;
        }
        const r = this[t];
        if (r === void 0) {
          console.warn(`THREE.Material: '${t}' is not a property of THREE.${this.type}.`);
          continue;
        }
        r && r.isColor ? r.set(n) : r && r.isVector3 && n && n.isVector3 ? r.copy(n) : this[t] = n;
      }
  }
  toJSON(e) {
    const t = e === void 0 || typeof e == "string";
    t && (e = {
      textures: {},
      images: {}
    });
    const n = {
      metadata: {
        version: 4.6,
        type: "Material",
        generator: "Material.toJSON"
      }
    };
    n.uuid = this.uuid, n.type = this.type, this.name !== "" && (n.name = this.name), this.color && this.color.isColor && (n.color = this.color.getHex()), this.roughness !== void 0 && (n.roughness = this.roughness), this.metalness !== void 0 && (n.metalness = this.metalness), this.sheen !== void 0 && (n.sheen = this.sheen), this.sheenColor && this.sheenColor.isColor && (n.sheenColor = this.sheenColor.getHex()), this.sheenRoughness !== void 0 && (n.sheenRoughness = this.sheenRoughness), this.emissive && this.emissive.isColor && (n.emissive = this.emissive.getHex()), this.emissiveIntensity && this.emissiveIntensity !== 1 && (n.emissiveIntensity = this.emissiveIntensity), this.specular && this.specular.isColor && (n.specular = this.specular.getHex()), this.specularIntensity !== void 0 && (n.specularIntensity = this.specularIntensity), this.specularColor && this.specularColor.isColor && (n.specularColor = this.specularColor.getHex()), this.shininess !== void 0 && (n.shininess = this.shininess), this.clearcoat !== void 0 && (n.clearcoat = this.clearcoat), this.clearcoatRoughness !== void 0 && (n.clearcoatRoughness = this.clearcoatRoughness), this.clearcoatMap && this.clearcoatMap.isTexture && (n.clearcoatMap = this.clearcoatMap.toJSON(e).uuid), this.clearcoatRoughnessMap && this.clearcoatRoughnessMap.isTexture && (n.clearcoatRoughnessMap = this.clearcoatRoughnessMap.toJSON(e).uuid), this.clearcoatNormalMap && this.clearcoatNormalMap.isTexture && (n.clearcoatNormalMap = this.clearcoatNormalMap.toJSON(e).uuid, n.clearcoatNormalScale = this.clearcoatNormalScale.toArray()), this.iridescence !== void 0 && (n.iridescence = this.iridescence), this.iridescenceIOR !== void 0 && (n.iridescenceIOR = this.iridescenceIOR), this.iridescenceThicknessRange !== void 0 && (n.iridescenceThicknessRange = this.iridescenceThicknessRange), this.iridescenceMap && this.iridescenceMap.isTexture && (n.iridescenceMap = this.iridescenceMap.toJSON(e).uuid), this.iridescenceThicknessMap && this.iridescenceThicknessMap.isTexture && (n.iridescenceThicknessMap = this.iridescenceThicknessMap.toJSON(e).uuid), this.anisotropy !== void 0 && (n.anisotropy = this.anisotropy), this.anisotropyRotation !== void 0 && (n.anisotropyRotation = this.anisotropyRotation), this.anisotropyMap && this.anisotropyMap.isTexture && (n.anisotropyMap = this.anisotropyMap.toJSON(e).uuid), this.map && this.map.isTexture && (n.map = this.map.toJSON(e).uuid), this.matcap && this.matcap.isTexture && (n.matcap = this.matcap.toJSON(e).uuid), this.alphaMap && this.alphaMap.isTexture && (n.alphaMap = this.alphaMap.toJSON(e).uuid), this.lightMap && this.lightMap.isTexture && (n.lightMap = this.lightMap.toJSON(e).uuid, n.lightMapIntensity = this.lightMapIntensity), this.aoMap && this.aoMap.isTexture && (n.aoMap = this.aoMap.toJSON(e).uuid, n.aoMapIntensity = this.aoMapIntensity), this.bumpMap && this.bumpMap.isTexture && (n.bumpMap = this.bumpMap.toJSON(e).uuid, n.bumpScale = this.bumpScale), this.normalMap && this.normalMap.isTexture && (n.normalMap = this.normalMap.toJSON(e).uuid, n.normalMapType = this.normalMapType, n.normalScale = this.normalScale.toArray()), this.displacementMap && this.displacementMap.isTexture && (n.displacementMap = this.displacementMap.toJSON(e).uuid, n.displacementScale = this.displacementScale, n.displacementBias = this.displacementBias), this.roughnessMap && this.roughnessMap.isTexture && (n.roughnessMap = this.roughnessMap.toJSON(e).uuid), this.metalnessMap && this.metalnessMap.isTexture && (n.metalnessMap = this.metalnessMap.toJSON(e).uuid), this.emissiveMap && this.emissiveMap.isTexture && (n.emissiveMap = this.emissiveMap.toJSON(e).uuid), this.specularMap && this.specularMap.isTexture && (n.specularMap = this.specularMap.toJSON(e).uuid), this.specularIntensityMap && this.specularIntensityMap.isTexture && (n.specularIntensityMap = this.specularIntensityMap.toJSON(e).uuid), this.specularColorMap && this.specularColorMap.isTexture && (n.specularColorMap = this.specularColorMap.toJSON(e).uuid), this.envMap && this.envMap.isTexture && (n.envMap = this.envMap.toJSON(e).uuid, this.combine !== void 0 && (n.combine = this.combine)), this.envMapIntensity !== void 0 && (n.envMapIntensity = this.envMapIntensity), this.reflectivity !== void 0 && (n.reflectivity = this.reflectivity), this.refractionRatio !== void 0 && (n.refractionRatio = this.refractionRatio), this.gradientMap && this.gradientMap.isTexture && (n.gradientMap = this.gradientMap.toJSON(e).uuid), this.transmission !== void 0 && (n.transmission = this.transmission), this.transmissionMap && this.transmissionMap.isTexture && (n.transmissionMap = this.transmissionMap.toJSON(e).uuid), this.thickness !== void 0 && (n.thickness = this.thickness), this.thicknessMap && this.thicknessMap.isTexture && (n.thicknessMap = this.thicknessMap.toJSON(e).uuid), this.attenuationDistance !== void 0 && this.attenuationDistance !== 1 / 0 && (n.attenuationDistance = this.attenuationDistance), this.attenuationColor !== void 0 && (n.attenuationColor = this.attenuationColor.getHex()), this.size !== void 0 && (n.size = this.size), this.shadowSide !== null && (n.shadowSide = this.shadowSide), this.sizeAttenuation !== void 0 && (n.sizeAttenuation = this.sizeAttenuation), this.blending !== $n && (n.blending = this.blending), this.side !== ln && (n.side = this.side), this.vertexColors && (n.vertexColors = !0), this.opacity < 1 && (n.opacity = this.opacity), this.transparent === !0 && (n.transparent = this.transparent), n.depthFunc = this.depthFunc, n.depthTest = this.depthTest, n.depthWrite = this.depthWrite, n.colorWrite = this.colorWrite, n.stencilWrite = this.stencilWrite, n.stencilWriteMask = this.stencilWriteMask, n.stencilFunc = this.stencilFunc, n.stencilRef = this.stencilRef, n.stencilFuncMask = this.stencilFuncMask, n.stencilFail = this.stencilFail, n.stencilZFail = this.stencilZFail, n.stencilZPass = this.stencilZPass, this.rotation !== void 0 && this.rotation !== 0 && (n.rotation = this.rotation), this.polygonOffset === !0 && (n.polygonOffset = !0), this.polygonOffsetFactor !== 0 && (n.polygonOffsetFactor = this.polygonOffsetFactor), this.polygonOffsetUnits !== 0 && (n.polygonOffsetUnits = this.polygonOffsetUnits), this.linewidth !== void 0 && this.linewidth !== 1 && (n.linewidth = this.linewidth), this.dashSize !== void 0 && (n.dashSize = this.dashSize), this.gapSize !== void 0 && (n.gapSize = this.gapSize), this.scale !== void 0 && (n.scale = this.scale), this.dithering === !0 && (n.dithering = !0), this.alphaTest > 0 && (n.alphaTest = this.alphaTest), this.alphaHash === !0 && (n.alphaHash = this.alphaHash), this.alphaToCoverage === !0 && (n.alphaToCoverage = this.alphaToCoverage), this.premultipliedAlpha === !0 && (n.premultipliedAlpha = this.premultipliedAlpha), this.forceSinglePass === !0 && (n.forceSinglePass = this.forceSinglePass), this.wireframe === !0 && (n.wireframe = this.wireframe), this.wireframeLinewidth > 1 && (n.wireframeLinewidth = this.wireframeLinewidth), this.wireframeLinecap !== "round" && (n.wireframeLinecap = this.wireframeLinecap), this.wireframeLinejoin !== "round" && (n.wireframeLinejoin = this.wireframeLinejoin), this.flatShading === !0 && (n.flatShading = this.flatShading), this.visible === !1 && (n.visible = !1), this.toneMapped === !1 && (n.toneMapped = !1), this.fog === !1 && (n.fog = !1), Object.keys(this.userData).length > 0 && (n.userData = this.userData);
    function r(s) {
      const o = [];
      for (const a in s) {
        const l = s[a];
        delete l.metadata, o.push(l);
      }
      return o;
    }
    if (t) {
      const s = r(e.textures), o = r(e.images);
      s.length > 0 && (n.textures = s), o.length > 0 && (n.images = o);
    }
    return n;
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    this.name = e.name, this.blending = e.blending, this.side = e.side, this.vertexColors = e.vertexColors, this.opacity = e.opacity, this.transparent = e.transparent, this.blendSrc = e.blendSrc, this.blendDst = e.blendDst, this.blendEquation = e.blendEquation, this.blendSrcAlpha = e.blendSrcAlpha, this.blendDstAlpha = e.blendDstAlpha, this.blendEquationAlpha = e.blendEquationAlpha, this.depthFunc = e.depthFunc, this.depthTest = e.depthTest, this.depthWrite = e.depthWrite, this.stencilWriteMask = e.stencilWriteMask, this.stencilFunc = e.stencilFunc, this.stencilRef = e.stencilRef, this.stencilFuncMask = e.stencilFuncMask, this.stencilFail = e.stencilFail, this.stencilZFail = e.stencilZFail, this.stencilZPass = e.stencilZPass, this.stencilWrite = e.stencilWrite;
    const t = e.clippingPlanes;
    let n = null;
    if (t !== null) {
      const r = t.length;
      n = new Array(r);
      for (let s = 0; s !== r; ++s)
        n[s] = t[s].clone();
    }
    return this.clippingPlanes = n, this.clipIntersection = e.clipIntersection, this.clipShadows = e.clipShadows, this.shadowSide = e.shadowSide, this.colorWrite = e.colorWrite, this.precision = e.precision, this.polygonOffset = e.polygonOffset, this.polygonOffsetFactor = e.polygonOffsetFactor, this.polygonOffsetUnits = e.polygonOffsetUnits, this.dithering = e.dithering, this.alphaTest = e.alphaTest, this.alphaHash = e.alphaHash, this.alphaToCoverage = e.alphaToCoverage, this.premultipliedAlpha = e.premultipliedAlpha, this.forceSinglePass = e.forceSinglePass, this.visible = e.visible, this.toneMapped = e.toneMapped, this.userData = JSON.parse(JSON.stringify(e.userData)), this;
  }
  dispose() {
    this.dispatchEvent({ type: "dispose" });
  }
  set needsUpdate(e) {
    e === !0 && this.version++;
  }
}
const _o = {
  aliceblue: 15792383,
  antiquewhite: 16444375,
  aqua: 65535,
  aquamarine: 8388564,
  azure: 15794175,
  beige: 16119260,
  bisque: 16770244,
  black: 0,
  blanchedalmond: 16772045,
  blue: 255,
  blueviolet: 9055202,
  brown: 10824234,
  burlywood: 14596231,
  cadetblue: 6266528,
  chartreuse: 8388352,
  chocolate: 13789470,
  coral: 16744272,
  cornflowerblue: 6591981,
  cornsilk: 16775388,
  crimson: 14423100,
  cyan: 65535,
  darkblue: 139,
  darkcyan: 35723,
  darkgoldenrod: 12092939,
  darkgray: 11119017,
  darkgreen: 25600,
  darkgrey: 11119017,
  darkkhaki: 12433259,
  darkmagenta: 9109643,
  darkolivegreen: 5597999,
  darkorange: 16747520,
  darkorchid: 10040012,
  darkred: 9109504,
  darksalmon: 15308410,
  darkseagreen: 9419919,
  darkslateblue: 4734347,
  darkslategray: 3100495,
  darkslategrey: 3100495,
  darkturquoise: 52945,
  darkviolet: 9699539,
  deeppink: 16716947,
  deepskyblue: 49151,
  dimgray: 6908265,
  dimgrey: 6908265,
  dodgerblue: 2003199,
  firebrick: 11674146,
  floralwhite: 16775920,
  forestgreen: 2263842,
  fuchsia: 16711935,
  gainsboro: 14474460,
  ghostwhite: 16316671,
  gold: 16766720,
  goldenrod: 14329120,
  gray: 8421504,
  green: 32768,
  greenyellow: 11403055,
  grey: 8421504,
  honeydew: 15794160,
  hotpink: 16738740,
  indianred: 13458524,
  indigo: 4915330,
  ivory: 16777200,
  khaki: 15787660,
  lavender: 15132410,
  lavenderblush: 16773365,
  lawngreen: 8190976,
  lemonchiffon: 16775885,
  lightblue: 11393254,
  lightcoral: 15761536,
  lightcyan: 14745599,
  lightgoldenrodyellow: 16448210,
  lightgray: 13882323,
  lightgreen: 9498256,
  lightgrey: 13882323,
  lightpink: 16758465,
  lightsalmon: 16752762,
  lightseagreen: 2142890,
  lightskyblue: 8900346,
  lightslategray: 7833753,
  lightslategrey: 7833753,
  lightsteelblue: 11584734,
  lightyellow: 16777184,
  lime: 65280,
  limegreen: 3329330,
  linen: 16445670,
  magenta: 16711935,
  maroon: 8388608,
  mediumaquamarine: 6737322,
  mediumblue: 205,
  mediumorchid: 12211667,
  mediumpurple: 9662683,
  mediumseagreen: 3978097,
  mediumslateblue: 8087790,
  mediumspringgreen: 64154,
  mediumturquoise: 4772300,
  mediumvioletred: 13047173,
  midnightblue: 1644912,
  mintcream: 16121850,
  mistyrose: 16770273,
  moccasin: 16770229,
  navajowhite: 16768685,
  navy: 128,
  oldlace: 16643558,
  olive: 8421376,
  olivedrab: 7048739,
  orange: 16753920,
  orangered: 16729344,
  orchid: 14315734,
  palegoldenrod: 15657130,
  palegreen: 10025880,
  paleturquoise: 11529966,
  palevioletred: 14381203,
  papayawhip: 16773077,
  peachpuff: 16767673,
  peru: 13468991,
  pink: 16761035,
  plum: 14524637,
  powderblue: 11591910,
  purple: 8388736,
  rebeccapurple: 6697881,
  red: 16711680,
  rosybrown: 12357519,
  royalblue: 4286945,
  saddlebrown: 9127187,
  salmon: 16416882,
  sandybrown: 16032864,
  seagreen: 3050327,
  seashell: 16774638,
  sienna: 10506797,
  silver: 12632256,
  skyblue: 8900331,
  slateblue: 6970061,
  slategray: 7372944,
  slategrey: 7372944,
  snow: 16775930,
  springgreen: 65407,
  steelblue: 4620980,
  tan: 13808780,
  teal: 32896,
  thistle: 14204888,
  tomato: 16737095,
  turquoise: 4251856,
  violet: 15631086,
  wheat: 16113331,
  white: 16777215,
  whitesmoke: 16119285,
  yellow: 16776960,
  yellowgreen: 10145074
}, Ct = { h: 0, s: 0, l: 0 }, Oi = { h: 0, s: 0, l: 0 };
function Dr(i, e, t) {
  return t < 0 && (t += 1), t > 1 && (t -= 1), t < 1 / 6 ? i + (e - i) * 6 * t : t < 1 / 2 ? e : t < 2 / 3 ? i + (e - i) * 6 * (2 / 3 - t) : i;
}
class We {
  constructor(e, t, n) {
    return this.isColor = !0, this.r = 1, this.g = 1, this.b = 1, this.set(e, t, n);
  }
  set(e, t, n) {
    if (t === void 0 && n === void 0) {
      const r = e;
      r && r.isColor ? this.copy(r) : typeof r == "number" ? this.setHex(r) : typeof r == "string" && this.setStyle(r);
    } else
      this.setRGB(e, t, n);
    return this;
  }
  setScalar(e) {
    return this.r = e, this.g = e, this.b = e, this;
  }
  setHex(e, t = Oe) {
    return e = Math.floor(e), this.r = (e >> 16 & 255) / 255, this.g = (e >> 8 & 255) / 255, this.b = (e & 255) / 255, bt.toWorkingColorSpace(this, t), this;
  }
  setRGB(e, t, n, r = bt.workingColorSpace) {
    return this.r = e, this.g = t, this.b = n, bt.toWorkingColorSpace(this, r), this;
  }
  setHSL(e, t, n, r = bt.workingColorSpace) {
    if (e = os(e, 1), t = it(t, 0, 1), n = it(n, 0, 1), t === 0)
      this.r = this.g = this.b = n;
    else {
      const s = n <= 0.5 ? n * (1 + t) : n + t - n * t, o = 2 * n - s;
      this.r = Dr(o, s, e + 1 / 3), this.g = Dr(o, s, e), this.b = Dr(o, s, e - 1 / 3);
    }
    return bt.toWorkingColorSpace(this, r), this;
  }
  setStyle(e, t = Oe) {
    function n(s) {
      s !== void 0 && parseFloat(s) < 1 && console.warn("THREE.Color: Alpha component of " + e + " will be ignored.");
    }
    let r;
    if (r = /^(\w+)\(([^\)]*)\)/.exec(e)) {
      let s;
      const o = r[1], a = r[2];
      switch (o) {
        case "rgb":
        case "rgba":
          if (s = /^\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(a))
            return n(s[4]), this.setRGB(
              Math.min(255, parseInt(s[1], 10)) / 255,
              Math.min(255, parseInt(s[2], 10)) / 255,
              Math.min(255, parseInt(s[3], 10)) / 255,
              t
            );
          if (s = /^\s*(\d+)\%\s*,\s*(\d+)\%\s*,\s*(\d+)\%\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(a))
            return n(s[4]), this.setRGB(
              Math.min(100, parseInt(s[1], 10)) / 100,
              Math.min(100, parseInt(s[2], 10)) / 100,
              Math.min(100, parseInt(s[3], 10)) / 100,
              t
            );
          break;
        case "hsl":
        case "hsla":
          if (s = /^\s*(\d*\.?\d+)\s*,\s*(\d*\.?\d+)\%\s*,\s*(\d*\.?\d+)\%\s*(?:,\s*(\d*\.?\d+)\s*)?$/.exec(a))
            return n(s[4]), this.setHSL(
              parseFloat(s[1]) / 360,
              parseFloat(s[2]) / 100,
              parseFloat(s[3]) / 100,
              t
            );
          break;
        default:
          console.warn("THREE.Color: Unknown color model " + e);
      }
    } else if (r = /^\#([A-Fa-f\d]+)$/.exec(e)) {
      const s = r[1], o = s.length;
      if (o === 3)
        return this.setRGB(
          parseInt(s.charAt(0), 16) / 15,
          parseInt(s.charAt(1), 16) / 15,
          parseInt(s.charAt(2), 16) / 15,
          t
        );
      if (o === 6)
        return this.setHex(parseInt(s, 16), t);
      console.warn("THREE.Color: Invalid hex color " + e);
    } else if (e && e.length > 0)
      return this.setColorName(e, t);
    return this;
  }
  setColorName(e, t = Oe) {
    const n = _o[e.toLowerCase()];
    return n !== void 0 ? this.setHex(n, t) : console.warn("THREE.Color: Unknown color " + e), this;
  }
  clone() {
    return new this.constructor(this.r, this.g, this.b);
  }
  copy(e) {
    return this.r = e.r, this.g = e.g, this.b = e.b, this;
  }
  copySRGBToLinear(e) {
    return this.r = Qn(e.r), this.g = Qn(e.g), this.b = Qn(e.b), this;
  }
  copyLinearToSRGB(e) {
    return this.r = Sr(e.r), this.g = Sr(e.g), this.b = Sr(e.b), this;
  }
  convertSRGBToLinear() {
    return this.copySRGBToLinear(this), this;
  }
  convertLinearToSRGB() {
    return this.copyLinearToSRGB(this), this;
  }
  getHex(e = Oe) {
    return bt.fromWorkingColorSpace(ct.copy(this), e), Math.round(it(ct.r * 255, 0, 255)) * 65536 + Math.round(it(ct.g * 255, 0, 255)) * 256 + Math.round(it(ct.b * 255, 0, 255));
  }
  getHexString(e = Oe) {
    return ("000000" + this.getHex(e).toString(16)).slice(-6);
  }
  getHSL(e, t = bt.workingColorSpace) {
    bt.fromWorkingColorSpace(ct.copy(this), t);
    const n = ct.r, r = ct.g, s = ct.b, o = Math.max(n, r, s), a = Math.min(n, r, s);
    let l, c;
    const h = (a + o) / 2;
    if (a === o)
      l = 0, c = 0;
    else {
      const f = o - a;
      switch (c = h <= 0.5 ? f / (o + a) : f / (2 - o - a), o) {
        case n:
          l = (r - s) / f + (r < s ? 6 : 0);
          break;
        case r:
          l = (s - n) / f + 2;
          break;
        case s:
          l = (n - r) / f + 4;
          break;
      }
      l /= 6;
    }
    return e.h = l, e.s = c, e.l = h, e;
  }
  getRGB(e, t = bt.workingColorSpace) {
    return bt.fromWorkingColorSpace(ct.copy(this), t), e.r = ct.r, e.g = ct.g, e.b = ct.b, e;
  }
  getStyle(e = Oe) {
    bt.fromWorkingColorSpace(ct.copy(this), e);
    const t = ct.r, n = ct.g, r = ct.b;
    return e !== Oe ? `color(${e} ${t.toFixed(3)} ${n.toFixed(3)} ${r.toFixed(3)})` : `rgb(${Math.round(t * 255)},${Math.round(n * 255)},${Math.round(r * 255)})`;
  }
  offsetHSL(e, t, n) {
    return this.getHSL(Ct), Ct.h += e, Ct.s += t, Ct.l += n, this.setHSL(Ct.h, Ct.s, Ct.l), this;
  }
  add(e) {
    return this.r += e.r, this.g += e.g, this.b += e.b, this;
  }
  addColors(e, t) {
    return this.r = e.r + t.r, this.g = e.g + t.g, this.b = e.b + t.b, this;
  }
  addScalar(e) {
    return this.r += e, this.g += e, this.b += e, this;
  }
  sub(e) {
    return this.r = Math.max(0, this.r - e.r), this.g = Math.max(0, this.g - e.g), this.b = Math.max(0, this.b - e.b), this;
  }
  multiply(e) {
    return this.r *= e.r, this.g *= e.g, this.b *= e.b, this;
  }
  multiplyScalar(e) {
    return this.r *= e, this.g *= e, this.b *= e, this;
  }
  lerp(e, t) {
    return this.r += (e.r - this.r) * t, this.g += (e.g - this.g) * t, this.b += (e.b - this.b) * t, this;
  }
  lerpColors(e, t, n) {
    return this.r = e.r + (t.r - e.r) * n, this.g = e.g + (t.g - e.g) * n, this.b = e.b + (t.b - e.b) * n, this;
  }
  lerpHSL(e, t) {
    this.getHSL(Ct), e.getHSL(Oi);
    const n = fi(Ct.h, Oi.h, t), r = fi(Ct.s, Oi.s, t), s = fi(Ct.l, Oi.l, t);
    return this.setHSL(n, r, s), this;
  }
  setFromVector3(e) {
    return this.r = e.x, this.g = e.y, this.b = e.z, this;
  }
  applyMatrix3(e) {
    const t = this.r, n = this.g, r = this.b, s = e.elements;
    return this.r = s[0] * t + s[3] * n + s[6] * r, this.g = s[1] * t + s[4] * n + s[7] * r, this.b = s[2] * t + s[5] * n + s[8] * r, this;
  }
  equals(e) {
    return e.r === this.r && e.g === this.g && e.b === this.b;
  }
  fromArray(e, t = 0) {
    return this.r = e[t], this.g = e[t + 1], this.b = e[t + 2], this;
  }
  toArray(e = [], t = 0) {
    return e[t] = this.r, e[t + 1] = this.g, e[t + 2] = this.b, e;
  }
  fromBufferAttribute(e, t) {
    return this.r = e.getX(t), this.g = e.getY(t), this.b = e.getZ(t), this;
  }
  toJSON() {
    return this.getHex();
  }
  *[Symbol.iterator]() {
    yield this.r, yield this.g, yield this.b;
  }
}
const ct = /* @__PURE__ */ new We();
We.NAMES = _o;
class vo extends Ai {
  constructor(e) {
    super(), this.isMeshBasicMaterial = !0, this.type = "MeshBasicMaterial", this.color = new We(16777215), this.map = null, this.lightMap = null, this.lightMapIntensity = 1, this.aoMap = null, this.aoMapIntensity = 1, this.specularMap = null, this.alphaMap = null, this.envMap = null, this.combine = to, this.reflectivity = 1, this.refractionRatio = 0.98, this.wireframe = !1, this.wireframeLinewidth = 1, this.wireframeLinecap = "round", this.wireframeLinejoin = "round", this.fog = !0, this.setValues(e);
  }
  copy(e) {
    return super.copy(e), this.color.copy(e.color), this.map = e.map, this.lightMap = e.lightMap, this.lightMapIntensity = e.lightMapIntensity, this.aoMap = e.aoMap, this.aoMapIntensity = e.aoMapIntensity, this.specularMap = e.specularMap, this.alphaMap = e.alphaMap, this.envMap = e.envMap, this.combine = e.combine, this.reflectivity = e.reflectivity, this.refractionRatio = e.refractionRatio, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.wireframeLinecap = e.wireframeLinecap, this.wireframeLinejoin = e.wireframeLinejoin, this.fog = e.fog, this;
  }
}
const tt = /* @__PURE__ */ new U(), Fi = /* @__PURE__ */ new oe();
class Ot {
  constructor(e, t, n = !1) {
    if (Array.isArray(e))
      throw new TypeError("THREE.BufferAttribute: array should be a Typed Array.");
    this.isBufferAttribute = !0, this.name = "", this.array = e, this.itemSize = t, this.count = e !== void 0 ? e.length / t : 0, this.normalized = n, this.usage = Qs, this.updateRange = { offset: 0, count: -1 }, this.gpuType = rn, this.version = 0;
  }
  onUploadCallback() {
  }
  set needsUpdate(e) {
    e === !0 && this.version++;
  }
  setUsage(e) {
    return this.usage = e, this;
  }
  copy(e) {
    return this.name = e.name, this.array = new e.array.constructor(e.array), this.itemSize = e.itemSize, this.count = e.count, this.normalized = e.normalized, this.usage = e.usage, this.gpuType = e.gpuType, this;
  }
  copyAt(e, t, n) {
    e *= this.itemSize, n *= t.itemSize;
    for (let r = 0, s = this.itemSize; r < s; r++)
      this.array[e + r] = t.array[n + r];
    return this;
  }
  copyArray(e) {
    return this.array.set(e), this;
  }
  applyMatrix3(e) {
    if (this.itemSize === 2)
      for (let t = 0, n = this.count; t < n; t++)
        Fi.fromBufferAttribute(this, t), Fi.applyMatrix3(e), this.setXY(t, Fi.x, Fi.y);
    else if (this.itemSize === 3)
      for (let t = 0, n = this.count; t < n; t++)
        tt.fromBufferAttribute(this, t), tt.applyMatrix3(e), this.setXYZ(t, tt.x, tt.y, tt.z);
    return this;
  }
  applyMatrix4(e) {
    for (let t = 0, n = this.count; t < n; t++)
      tt.fromBufferAttribute(this, t), tt.applyMatrix4(e), this.setXYZ(t, tt.x, tt.y, tt.z);
    return this;
  }
  applyNormalMatrix(e) {
    for (let t = 0, n = this.count; t < n; t++)
      tt.fromBufferAttribute(this, t), tt.applyNormalMatrix(e), this.setXYZ(t, tt.x, tt.y, tt.z);
    return this;
  }
  transformDirection(e) {
    for (let t = 0, n = this.count; t < n; t++)
      tt.fromBufferAttribute(this, t), tt.transformDirection(e), this.setXYZ(t, tt.x, tt.y, tt.z);
    return this;
  }
  set(e, t = 0) {
    return this.array.set(e, t), this;
  }
  getComponent(e, t) {
    let n = this.array[e * this.itemSize + t];
    return this.normalized && (n = Zn(n, this.array)), n;
  }
  setComponent(e, t, n) {
    return this.normalized && (n = dt(n, this.array)), this.array[e * this.itemSize + t] = n, this;
  }
  getX(e) {
    let t = this.array[e * this.itemSize];
    return this.normalized && (t = Zn(t, this.array)), t;
  }
  setX(e, t) {
    return this.normalized && (t = dt(t, this.array)), this.array[e * this.itemSize] = t, this;
  }
  getY(e) {
    let t = this.array[e * this.itemSize + 1];
    return this.normalized && (t = Zn(t, this.array)), t;
  }
  setY(e, t) {
    return this.normalized && (t = dt(t, this.array)), this.array[e * this.itemSize + 1] = t, this;
  }
  getZ(e) {
    let t = this.array[e * this.itemSize + 2];
    return this.normalized && (t = Zn(t, this.array)), t;
  }
  setZ(e, t) {
    return this.normalized && (t = dt(t, this.array)), this.array[e * this.itemSize + 2] = t, this;
  }
  getW(e) {
    let t = this.array[e * this.itemSize + 3];
    return this.normalized && (t = Zn(t, this.array)), t;
  }
  setW(e, t) {
    return this.normalized && (t = dt(t, this.array)), this.array[e * this.itemSize + 3] = t, this;
  }
  setXY(e, t, n) {
    return e *= this.itemSize, this.normalized && (t = dt(t, this.array), n = dt(n, this.array)), this.array[e + 0] = t, this.array[e + 1] = n, this;
  }
  setXYZ(e, t, n, r) {
    return e *= this.itemSize, this.normalized && (t = dt(t, this.array), n = dt(n, this.array), r = dt(r, this.array)), this.array[e + 0] = t, this.array[e + 1] = n, this.array[e + 2] = r, this;
  }
  setXYZW(e, t, n, r, s) {
    return e *= this.itemSize, this.normalized && (t = dt(t, this.array), n = dt(n, this.array), r = dt(r, this.array), s = dt(s, this.array)), this.array[e + 0] = t, this.array[e + 1] = n, this.array[e + 2] = r, this.array[e + 3] = s, this;
  }
  onUpload(e) {
    return this.onUploadCallback = e, this;
  }
  clone() {
    return new this.constructor(this.array, this.itemSize).copy(this);
  }
  toJSON() {
    const e = {
      itemSize: this.itemSize,
      type: this.array.constructor.name,
      array: Array.from(this.array),
      normalized: this.normalized
    };
    return this.name !== "" && (e.name = this.name), this.usage !== Qs && (e.usage = this.usage), (this.updateRange.offset !== 0 || this.updateRange.count !== -1) && (e.updateRange = this.updateRange), e;
  }
}
class xo extends Ot {
  constructor(e, t, n) {
    super(new Uint16Array(e), t, n);
  }
}
class Mo extends Ot {
  constructor(e, t, n) {
    super(new Uint32Array(e), t, n);
  }
}
class jt extends Ot {
  constructor(e, t, n) {
    super(new Float32Array(e), t, n);
  }
}
let Ac = 0;
const yt = /* @__PURE__ */ new nt(), Ir = /* @__PURE__ */ new ht(), Gn = /* @__PURE__ */ new U(), xt = /* @__PURE__ */ new Ti(), ci = /* @__PURE__ */ new Ti(), at = /* @__PURE__ */ new U();
class cn extends Rn {
  constructor() {
    super(), this.isBufferGeometry = !0, Object.defineProperty(this, "id", { value: Ac++ }), this.uuid = Cn(), this.name = "", this.type = "BufferGeometry", this.index = null, this.attributes = {}, this.morphAttributes = {}, this.morphTargetsRelative = !1, this.groups = [], this.boundingBox = null, this.boundingSphere = null, this.drawRange = { start: 0, count: 1 / 0 }, this.userData = {};
  }
  getIndex() {
    return this.index;
  }
  setIndex(e) {
    return Array.isArray(e) ? this.index = new (fo(e) ? Mo : xo)(e, 1) : this.index = e, this;
  }
  getAttribute(e) {
    return this.attributes[e];
  }
  setAttribute(e, t) {
    return this.attributes[e] = t, this;
  }
  deleteAttribute(e) {
    return delete this.attributes[e], this;
  }
  hasAttribute(e) {
    return this.attributes[e] !== void 0;
  }
  addGroup(e, t, n = 0) {
    this.groups.push({
      start: e,
      count: t,
      materialIndex: n
    });
  }
  clearGroups() {
    this.groups = [];
  }
  setDrawRange(e, t) {
    this.drawRange.start = e, this.drawRange.count = t;
  }
  applyMatrix4(e) {
    const t = this.attributes.position;
    t !== void 0 && (t.applyMatrix4(e), t.needsUpdate = !0);
    const n = this.attributes.normal;
    if (n !== void 0) {
      const s = new Be().getNormalMatrix(e);
      n.applyNormalMatrix(s), n.needsUpdate = !0;
    }
    const r = this.attributes.tangent;
    return r !== void 0 && (r.transformDirection(e), r.needsUpdate = !0), this.boundingBox !== null && this.computeBoundingBox(), this.boundingSphere !== null && this.computeBoundingSphere(), this;
  }
  applyQuaternion(e) {
    return yt.makeRotationFromQuaternion(e), this.applyMatrix4(yt), this;
  }
  rotateX(e) {
    return yt.makeRotationX(e), this.applyMatrix4(yt), this;
  }
  rotateY(e) {
    return yt.makeRotationY(e), this.applyMatrix4(yt), this;
  }
  rotateZ(e) {
    return yt.makeRotationZ(e), this.applyMatrix4(yt), this;
  }
  translate(e, t, n) {
    return yt.makeTranslation(e, t, n), this.applyMatrix4(yt), this;
  }
  scale(e, t, n) {
    return yt.makeScale(e, t, n), this.applyMatrix4(yt), this;
  }
  lookAt(e) {
    return Ir.lookAt(e), Ir.updateMatrix(), this.applyMatrix4(Ir.matrix), this;
  }
  center() {
    return this.computeBoundingBox(), this.boundingBox.getCenter(Gn).negate(), this.translate(Gn.x, Gn.y, Gn.z), this;
  }
  setFromPoints(e) {
    const t = [];
    for (let n = 0, r = e.length; n < r; n++) {
      const s = e[n];
      t.push(s.x, s.y, s.z || 0);
    }
    return this.setAttribute("position", new jt(t, 3)), this;
  }
  computeBoundingBox() {
    this.boundingBox === null && (this.boundingBox = new Ti());
    const e = this.attributes.position, t = this.morphAttributes.position;
    if (e && e.isGLBufferAttribute) {
      console.error('THREE.BufferGeometry.computeBoundingBox(): GLBufferAttribute requires a manual bounding box. Alternatively set "mesh.frustumCulled" to "false".', this), this.boundingBox.set(
        new U(-1 / 0, -1 / 0, -1 / 0),
        new U(1 / 0, 1 / 0, 1 / 0)
      );
      return;
    }
    if (e !== void 0) {
      if (this.boundingBox.setFromBufferAttribute(e), t)
        for (let n = 0, r = t.length; n < r; n++) {
          const s = t[n];
          xt.setFromBufferAttribute(s), this.morphTargetsRelative ? (at.addVectors(this.boundingBox.min, xt.min), this.boundingBox.expandByPoint(at), at.addVectors(this.boundingBox.max, xt.max), this.boundingBox.expandByPoint(at)) : (this.boundingBox.expandByPoint(xt.min), this.boundingBox.expandByPoint(xt.max));
        }
    } else
      this.boundingBox.makeEmpty();
    (isNaN(this.boundingBox.min.x) || isNaN(this.boundingBox.min.y) || isNaN(this.boundingBox.min.z)) && console.error('THREE.BufferGeometry.computeBoundingBox(): Computed min/max have NaN values. The "position" attribute is likely to have NaN values.', this);
  }
  computeBoundingSphere() {
    this.boundingSphere === null && (this.boundingSphere = new ls());
    const e = this.attributes.position, t = this.morphAttributes.position;
    if (e && e.isGLBufferAttribute) {
      console.error('THREE.BufferGeometry.computeBoundingSphere(): GLBufferAttribute requires a manual bounding sphere. Alternatively set "mesh.frustumCulled" to "false".', this), this.boundingSphere.set(new U(), 1 / 0);
      return;
    }
    if (e) {
      const n = this.boundingSphere.center;
      if (xt.setFromBufferAttribute(e), t)
        for (let s = 0, o = t.length; s < o; s++) {
          const a = t[s];
          ci.setFromBufferAttribute(a), this.morphTargetsRelative ? (at.addVectors(xt.min, ci.min), xt.expandByPoint(at), at.addVectors(xt.max, ci.max), xt.expandByPoint(at)) : (xt.expandByPoint(ci.min), xt.expandByPoint(ci.max));
        }
      xt.getCenter(n);
      let r = 0;
      for (let s = 0, o = e.count; s < o; s++)
        at.fromBufferAttribute(e, s), r = Math.max(r, n.distanceToSquared(at));
      if (t)
        for (let s = 0, o = t.length; s < o; s++) {
          const a = t[s], l = this.morphTargetsRelative;
          for (let c = 0, h = a.count; c < h; c++)
            at.fromBufferAttribute(a, c), l && (Gn.fromBufferAttribute(e, c), at.add(Gn)), r = Math.max(r, n.distanceToSquared(at));
        }
      this.boundingSphere.radius = Math.sqrt(r), isNaN(this.boundingSphere.radius) && console.error('THREE.BufferGeometry.computeBoundingSphere(): Computed radius is NaN. The "position" attribute is likely to have NaN values.', this);
    }
  }
  computeTangents() {
    const e = this.index, t = this.attributes;
    if (e === null || t.position === void 0 || t.normal === void 0 || t.uv === void 0) {
      console.error("THREE.BufferGeometry: .computeTangents() failed. Missing required attributes (index, position, normal or uv)");
      return;
    }
    const n = e.array, r = t.position.array, s = t.normal.array, o = t.uv.array, a = r.length / 3;
    this.hasAttribute("tangent") === !1 && this.setAttribute("tangent", new Ot(new Float32Array(4 * a), 4));
    const l = this.getAttribute("tangent").array, c = [], h = [];
    for (let y = 0; y < a; y++)
      c[y] = new U(), h[y] = new U();
    const f = new U(), u = new U(), m = new U(), g = new oe(), x = new oe(), p = new oe(), d = new U(), A = new U();
    function _(y, Y, ce) {
      f.fromArray(r, y * 3), u.fromArray(r, Y * 3), m.fromArray(r, ce * 3), g.fromArray(o, y * 2), x.fromArray(o, Y * 2), p.fromArray(o, ce * 2), u.sub(f), m.sub(f), x.sub(g), p.sub(g);
      const B = 1 / (x.x * p.y - p.x * x.y);
      isFinite(B) && (d.copy(u).multiplyScalar(p.y).addScaledVector(m, -x.y).multiplyScalar(B), A.copy(m).multiplyScalar(x.x).addScaledVector(u, -p.x).multiplyScalar(B), c[y].add(d), c[Y].add(d), c[ce].add(d), h[y].add(A), h[Y].add(A), h[ce].add(A));
    }
    let T = this.groups;
    T.length === 0 && (T = [{
      start: 0,
      count: n.length
    }]);
    for (let y = 0, Y = T.length; y < Y; ++y) {
      const ce = T[y], B = ce.start, H = ce.count;
      for (let G = B, Q = B + H; G < Q; G += 3)
        _(
          n[G + 0],
          n[G + 1],
          n[G + 2]
        );
    }
    const C = new U(), L = new U(), w = new U(), V = new U();
    function M(y) {
      w.fromArray(s, y * 3), V.copy(w);
      const Y = c[y];
      C.copy(Y), C.sub(w.multiplyScalar(w.dot(Y))).normalize(), L.crossVectors(V, Y);
      const B = L.dot(h[y]) < 0 ? -1 : 1;
      l[y * 4] = C.x, l[y * 4 + 1] = C.y, l[y * 4 + 2] = C.z, l[y * 4 + 3] = B;
    }
    for (let y = 0, Y = T.length; y < Y; ++y) {
      const ce = T[y], B = ce.start, H = ce.count;
      for (let G = B, Q = B + H; G < Q; G += 3)
        M(n[G + 0]), M(n[G + 1]), M(n[G + 2]);
    }
  }
  computeVertexNormals() {
    const e = this.index, t = this.getAttribute("position");
    if (t !== void 0) {
      let n = this.getAttribute("normal");
      if (n === void 0)
        n = new Ot(new Float32Array(t.count * 3), 3), this.setAttribute("normal", n);
      else
        for (let u = 0, m = n.count; u < m; u++)
          n.setXYZ(u, 0, 0, 0);
      const r = new U(), s = new U(), o = new U(), a = new U(), l = new U(), c = new U(), h = new U(), f = new U();
      if (e)
        for (let u = 0, m = e.count; u < m; u += 3) {
          const g = e.getX(u + 0), x = e.getX(u + 1), p = e.getX(u + 2);
          r.fromBufferAttribute(t, g), s.fromBufferAttribute(t, x), o.fromBufferAttribute(t, p), h.subVectors(o, s), f.subVectors(r, s), h.cross(f), a.fromBufferAttribute(n, g), l.fromBufferAttribute(n, x), c.fromBufferAttribute(n, p), a.add(h), l.add(h), c.add(h), n.setXYZ(g, a.x, a.y, a.z), n.setXYZ(x, l.x, l.y, l.z), n.setXYZ(p, c.x, c.y, c.z);
        }
      else
        for (let u = 0, m = t.count; u < m; u += 3)
          r.fromBufferAttribute(t, u + 0), s.fromBufferAttribute(t, u + 1), o.fromBufferAttribute(t, u + 2), h.subVectors(o, s), f.subVectors(r, s), h.cross(f), n.setXYZ(u + 0, h.x, h.y, h.z), n.setXYZ(u + 1, h.x, h.y, h.z), n.setXYZ(u + 2, h.x, h.y, h.z);
      this.normalizeNormals(), n.needsUpdate = !0;
    }
  }
  normalizeNormals() {
    const e = this.attributes.normal;
    for (let t = 0, n = e.count; t < n; t++)
      at.fromBufferAttribute(e, t), at.normalize(), e.setXYZ(t, at.x, at.y, at.z);
  }
  toNonIndexed() {
    function e(a, l) {
      const c = a.array, h = a.itemSize, f = a.normalized, u = new c.constructor(l.length * h);
      let m = 0, g = 0;
      for (let x = 0, p = l.length; x < p; x++) {
        a.isInterleavedBufferAttribute ? m = l[x] * a.data.stride + a.offset : m = l[x] * h;
        for (let d = 0; d < h; d++)
          u[g++] = c[m++];
      }
      return new Ot(u, h, f);
    }
    if (this.index === null)
      return console.warn("THREE.BufferGeometry.toNonIndexed(): BufferGeometry is already non-indexed."), this;
    const t = new cn(), n = this.index.array, r = this.attributes;
    for (const a in r) {
      const l = r[a], c = e(l, n);
      t.setAttribute(a, c);
    }
    const s = this.morphAttributes;
    for (const a in s) {
      const l = [], c = s[a];
      for (let h = 0, f = c.length; h < f; h++) {
        const u = c[h], m = e(u, n);
        l.push(m);
      }
      t.morphAttributes[a] = l;
    }
    t.morphTargetsRelative = this.morphTargetsRelative;
    const o = this.groups;
    for (let a = 0, l = o.length; a < l; a++) {
      const c = o[a];
      t.addGroup(c.start, c.count, c.materialIndex);
    }
    return t;
  }
  toJSON() {
    const e = {
      metadata: {
        version: 4.6,
        type: "BufferGeometry",
        generator: "BufferGeometry.toJSON"
      }
    };
    if (e.uuid = this.uuid, e.type = this.type, this.name !== "" && (e.name = this.name), Object.keys(this.userData).length > 0 && (e.userData = this.userData), this.parameters !== void 0) {
      const l = this.parameters;
      for (const c in l)
        l[c] !== void 0 && (e[c] = l[c]);
      return e;
    }
    e.data = { attributes: {} };
    const t = this.index;
    t !== null && (e.data.index = {
      type: t.array.constructor.name,
      array: Array.prototype.slice.call(t.array)
    });
    const n = this.attributes;
    for (const l in n) {
      const c = n[l];
      e.data.attributes[l] = c.toJSON(e.data);
    }
    const r = {};
    let s = !1;
    for (const l in this.morphAttributes) {
      const c = this.morphAttributes[l], h = [];
      for (let f = 0, u = c.length; f < u; f++) {
        const m = c[f];
        h.push(m.toJSON(e.data));
      }
      h.length > 0 && (r[l] = h, s = !0);
    }
    s && (e.data.morphAttributes = r, e.data.morphTargetsRelative = this.morphTargetsRelative);
    const o = this.groups;
    o.length > 0 && (e.data.groups = JSON.parse(JSON.stringify(o)));
    const a = this.boundingSphere;
    return a !== null && (e.data.boundingSphere = {
      center: a.center.toArray(),
      radius: a.radius
    }), e;
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    this.index = null, this.attributes = {}, this.morphAttributes = {}, this.groups = [], this.boundingBox = null, this.boundingSphere = null;
    const t = {};
    this.name = e.name;
    const n = e.index;
    n !== null && this.setIndex(n.clone(t));
    const r = e.attributes;
    for (const c in r) {
      const h = r[c];
      this.setAttribute(c, h.clone(t));
    }
    const s = e.morphAttributes;
    for (const c in s) {
      const h = [], f = s[c];
      for (let u = 0, m = f.length; u < m; u++)
        h.push(f[u].clone(t));
      this.morphAttributes[c] = h;
    }
    this.morphTargetsRelative = e.morphTargetsRelative;
    const o = e.groups;
    for (let c = 0, h = o.length; c < h; c++) {
      const f = o[c];
      this.addGroup(f.start, f.count, f.materialIndex);
    }
    const a = e.boundingBox;
    a !== null && (this.boundingBox = a.clone());
    const l = e.boundingSphere;
    return l !== null && (this.boundingSphere = l.clone()), this.drawRange.start = e.drawRange.start, this.drawRange.count = e.drawRange.count, this.userData = e.userData, this;
  }
  dispose() {
    this.dispatchEvent({ type: "dispose" });
  }
}
const fa = /* @__PURE__ */ new nt(), pn = /* @__PURE__ */ new cs(), Bi = /* @__PURE__ */ new ls(), da = /* @__PURE__ */ new U(), Vn = /* @__PURE__ */ new U(), kn = /* @__PURE__ */ new U(), Wn = /* @__PURE__ */ new U(), Nr = /* @__PURE__ */ new U(), zi = /* @__PURE__ */ new U(), Hi = /* @__PURE__ */ new oe(), Gi = /* @__PURE__ */ new oe(), Vi = /* @__PURE__ */ new oe(), pa = /* @__PURE__ */ new U(), ma = /* @__PURE__ */ new U(), ga = /* @__PURE__ */ new U(), ki = /* @__PURE__ */ new U(), Wi = /* @__PURE__ */ new U();
class Nt extends ht {
  constructor(e = new cn(), t = new vo()) {
    super(), this.isMesh = !0, this.type = "Mesh", this.geometry = e, this.material = t, this.updateMorphTargets();
  }
  copy(e, t) {
    return super.copy(e, t), e.morphTargetInfluences !== void 0 && (this.morphTargetInfluences = e.morphTargetInfluences.slice()), e.morphTargetDictionary !== void 0 && (this.morphTargetDictionary = Object.assign({}, e.morphTargetDictionary)), this.material = e.material, this.geometry = e.geometry, this;
  }
  updateMorphTargets() {
    const t = this.geometry.morphAttributes, n = Object.keys(t);
    if (n.length > 0) {
      const r = t[n[0]];
      if (r !== void 0) {
        this.morphTargetInfluences = [], this.morphTargetDictionary = {};
        for (let s = 0, o = r.length; s < o; s++) {
          const a = r[s].name || String(s);
          this.morphTargetInfluences.push(0), this.morphTargetDictionary[a] = s;
        }
      }
    }
  }
  getVertexPosition(e, t) {
    const n = this.geometry, r = n.attributes.position, s = n.morphAttributes.position, o = n.morphTargetsRelative;
    t.fromBufferAttribute(r, e);
    const a = this.morphTargetInfluences;
    if (s && a) {
      zi.set(0, 0, 0);
      for (let l = 0, c = s.length; l < c; l++) {
        const h = a[l], f = s[l];
        h !== 0 && (Nr.fromBufferAttribute(f, e), o ? zi.addScaledVector(Nr, h) : zi.addScaledVector(Nr.sub(t), h));
      }
      t.add(zi);
    }
    return t;
  }
  raycast(e, t) {
    const n = this.geometry, r = this.material, s = this.matrixWorld;
    r !== void 0 && (n.boundingSphere === null && n.computeBoundingSphere(), Bi.copy(n.boundingSphere), Bi.applyMatrix4(s), pn.copy(e.ray).recast(e.near), !(Bi.containsPoint(pn.origin) === !1 && (pn.intersectSphere(Bi, da) === null || pn.origin.distanceToSquared(da) > (e.far - e.near) ** 2)) && (fa.copy(s).invert(), pn.copy(e.ray).applyMatrix4(fa), !(n.boundingBox !== null && pn.intersectsBox(n.boundingBox) === !1) && this._computeIntersections(e, t, pn)));
  }
  _computeIntersections(e, t, n) {
    let r;
    const s = this.geometry, o = this.material, a = s.index, l = s.attributes.position, c = s.attributes.uv, h = s.attributes.uv1, f = s.attributes.normal, u = s.groups, m = s.drawRange;
    if (a !== null)
      if (Array.isArray(o))
        for (let g = 0, x = u.length; g < x; g++) {
          const p = u[g], d = o[p.materialIndex], A = Math.max(p.start, m.start), _ = Math.min(a.count, Math.min(p.start + p.count, m.start + m.count));
          for (let T = A, C = _; T < C; T += 3) {
            const L = a.getX(T), w = a.getX(T + 1), V = a.getX(T + 2);
            r = Xi(this, d, e, n, c, h, f, L, w, V), r && (r.faceIndex = Math.floor(T / 3), r.face.materialIndex = p.materialIndex, t.push(r));
          }
        }
      else {
        const g = Math.max(0, m.start), x = Math.min(a.count, m.start + m.count);
        for (let p = g, d = x; p < d; p += 3) {
          const A = a.getX(p), _ = a.getX(p + 1), T = a.getX(p + 2);
          r = Xi(this, o, e, n, c, h, f, A, _, T), r && (r.faceIndex = Math.floor(p / 3), t.push(r));
        }
      }
    else if (l !== void 0)
      if (Array.isArray(o))
        for (let g = 0, x = u.length; g < x; g++) {
          const p = u[g], d = o[p.materialIndex], A = Math.max(p.start, m.start), _ = Math.min(l.count, Math.min(p.start + p.count, m.start + m.count));
          for (let T = A, C = _; T < C; T += 3) {
            const L = T, w = T + 1, V = T + 2;
            r = Xi(this, d, e, n, c, h, f, L, w, V), r && (r.faceIndex = Math.floor(T / 3), r.face.materialIndex = p.materialIndex, t.push(r));
          }
        }
      else {
        const g = Math.max(0, m.start), x = Math.min(l.count, m.start + m.count);
        for (let p = g, d = x; p < d; p += 3) {
          const A = p, _ = p + 1, T = p + 2;
          r = Xi(this, o, e, n, c, h, f, A, _, T), r && (r.faceIndex = Math.floor(p / 3), t.push(r));
        }
      }
  }
}
function bc(i, e, t, n, r, s, o, a) {
  let l;
  if (e.side === gt ? l = n.intersectTriangle(o, s, r, !0, a) : l = n.intersectTriangle(r, s, o, e.side === ln, a), l === null)
    return null;
  Wi.copy(a), Wi.applyMatrix4(i.matrixWorld);
  const c = t.ray.origin.distanceTo(Wi);
  return c < t.near || c > t.far ? null : {
    distance: c,
    point: Wi.clone(),
    object: i
  };
}
function Xi(i, e, t, n, r, s, o, a, l, c) {
  i.getVertexPosition(a, Vn), i.getVertexPosition(l, kn), i.getVertexPosition(c, Wn);
  const h = bc(i, e, t, n, Vn, kn, Wn, ki);
  if (h) {
    r && (Hi.fromBufferAttribute(r, a), Gi.fromBufferAttribute(r, l), Vi.fromBufferAttribute(r, c), h.uv = Pt.getInterpolation(ki, Vn, kn, Wn, Hi, Gi, Vi, new oe())), s && (Hi.fromBufferAttribute(s, a), Gi.fromBufferAttribute(s, l), Vi.fromBufferAttribute(s, c), h.uv1 = Pt.getInterpolation(ki, Vn, kn, Wn, Hi, Gi, Vi, new oe()), h.uv2 = h.uv1), o && (pa.fromBufferAttribute(o, a), ma.fromBufferAttribute(o, l), ga.fromBufferAttribute(o, c), h.normal = Pt.getInterpolation(ki, Vn, kn, Wn, pa, ma, ga, new U()), h.normal.dot(n.direction) > 0 && h.normal.multiplyScalar(-1));
    const f = {
      a,
      b: l,
      c,
      normal: new U(),
      materialIndex: 0
    };
    Pt.getNormal(Vn, kn, Wn, f.normal), h.face = f;
  }
  return h;
}
class bi extends cn {
  constructor(e = 1, t = 1, n = 1, r = 1, s = 1, o = 1) {
    super(), this.type = "BoxGeometry", this.parameters = {
      width: e,
      height: t,
      depth: n,
      widthSegments: r,
      heightSegments: s,
      depthSegments: o
    };
    const a = this;
    r = Math.floor(r), s = Math.floor(s), o = Math.floor(o);
    const l = [], c = [], h = [], f = [];
    let u = 0, m = 0;
    g("z", "y", "x", -1, -1, n, t, e, o, s, 0), g("z", "y", "x", 1, -1, n, t, -e, o, s, 1), g("x", "z", "y", 1, 1, e, n, t, r, o, 2), g("x", "z", "y", 1, -1, e, n, -t, r, o, 3), g("x", "y", "z", 1, -1, e, t, n, r, s, 4), g("x", "y", "z", -1, -1, e, t, -n, r, s, 5), this.setIndex(l), this.setAttribute("position", new jt(c, 3)), this.setAttribute("normal", new jt(h, 3)), this.setAttribute("uv", new jt(f, 2));
    function g(x, p, d, A, _, T, C, L, w, V, M) {
      const y = T / w, Y = C / V, ce = T / 2, B = C / 2, H = L / 2, G = w + 1, Q = V + 1;
      let X = 0, j = 0;
      const J = new U();
      for (let ee = 0; ee < Q; ee++) {
        const I = ee * Y - B;
        for (let q = 0; q < G; q++) {
          const pe = q * y - ce;
          J[x] = pe * A, J[p] = I * _, J[d] = H, c.push(J.x, J.y, J.z), J[x] = 0, J[p] = 0, J[d] = L > 0 ? 1 : -1, h.push(J.x, J.y, J.z), f.push(q / w), f.push(1 - ee / V), X += 1;
        }
      }
      for (let ee = 0; ee < V; ee++)
        for (let I = 0; I < w; I++) {
          const q = u + I + G * ee, pe = u + I + G * (ee + 1), me = u + (I + 1) + G * (ee + 1), ve = u + (I + 1) + G * ee;
          l.push(q, pe, ve), l.push(pe, me, ve), j += 6;
        }
      a.addGroup(m, j, M), m += j, u += X;
    }
  }
  copy(e) {
    return super.copy(e), this.parameters = Object.assign({}, e.parameters), this;
  }
  static fromJSON(e) {
    return new bi(e.width, e.height, e.depth, e.widthSegments, e.heightSegments, e.depthSegments);
  }
}
function ii(i) {
  const e = {};
  for (const t in i) {
    e[t] = {};
    for (const n in i[t]) {
      const r = i[t][n];
      r && (r.isColor || r.isMatrix3 || r.isMatrix4 || r.isVector2 || r.isVector3 || r.isVector4 || r.isTexture || r.isQuaternion) ? r.isRenderTargetTexture ? (console.warn("UniformsUtils: Textures of render targets cannot be cloned via cloneUniforms() or mergeUniforms()."), e[t][n] = null) : e[t][n] = r.clone() : Array.isArray(r) ? e[t][n] = r.slice() : e[t][n] = r;
    }
  }
  return e;
}
function pt(i) {
  const e = {};
  for (let t = 0; t < i.length; t++) {
    const n = ii(i[t]);
    for (const r in n)
      e[r] = n[r];
  }
  return e;
}
function wc(i) {
  const e = [];
  for (let t = 0; t < i.length; t++)
    e.push(i[t].clone());
  return e;
}
function So(i) {
  return i.getRenderTarget() === null ? i.outputColorSpace : Ft;
}
const Rc = { clone: ii, merge: pt };
var Cc = `void main() {
	gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );
}`, Pc = `void main() {
	gl_FragColor = vec4( 1.0, 0.0, 0.0, 1.0 );
}`;
class bn extends Ai {
  constructor(e) {
    super(), this.isShaderMaterial = !0, this.type = "ShaderMaterial", this.defines = {}, this.uniforms = {}, this.uniformsGroups = [], this.vertexShader = Cc, this.fragmentShader = Pc, this.linewidth = 1, this.wireframe = !1, this.wireframeLinewidth = 1, this.fog = !1, this.lights = !1, this.clipping = !1, this.forceSinglePass = !0, this.extensions = {
      derivatives: !1,
      // set to use derivatives
      fragDepth: !1,
      // set to use fragment depth values
      drawBuffers: !1,
      // set to use draw buffers
      shaderTextureLOD: !1
      // set to use shader texture LOD
    }, this.defaultAttributeValues = {
      color: [1, 1, 1],
      uv: [0, 0],
      uv1: [0, 0]
    }, this.index0AttributeName = void 0, this.uniformsNeedUpdate = !1, this.glslVersion = null, e !== void 0 && this.setValues(e);
  }
  copy(e) {
    return super.copy(e), this.fragmentShader = e.fragmentShader, this.vertexShader = e.vertexShader, this.uniforms = ii(e.uniforms), this.uniformsGroups = wc(e.uniformsGroups), this.defines = Object.assign({}, e.defines), this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.fog = e.fog, this.lights = e.lights, this.clipping = e.clipping, this.extensions = Object.assign({}, e.extensions), this.glslVersion = e.glslVersion, this;
  }
  toJSON(e) {
    const t = super.toJSON(e);
    t.glslVersion = this.glslVersion, t.uniforms = {};
    for (const r in this.uniforms) {
      const o = this.uniforms[r].value;
      o && o.isTexture ? t.uniforms[r] = {
        type: "t",
        value: o.toJSON(e).uuid
      } : o && o.isColor ? t.uniforms[r] = {
        type: "c",
        value: o.getHex()
      } : o && o.isVector2 ? t.uniforms[r] = {
        type: "v2",
        value: o.toArray()
      } : o && o.isVector3 ? t.uniforms[r] = {
        type: "v3",
        value: o.toArray()
      } : o && o.isVector4 ? t.uniforms[r] = {
        type: "v4",
        value: o.toArray()
      } : o && o.isMatrix3 ? t.uniforms[r] = {
        type: "m3",
        value: o.toArray()
      } : o && o.isMatrix4 ? t.uniforms[r] = {
        type: "m4",
        value: o.toArray()
      } : t.uniforms[r] = {
        value: o
      };
    }
    Object.keys(this.defines).length > 0 && (t.defines = this.defines), t.vertexShader = this.vertexShader, t.fragmentShader = this.fragmentShader, t.lights = this.lights, t.clipping = this.clipping;
    const n = {};
    for (const r in this.extensions)
      this.extensions[r] === !0 && (n[r] = !0);
    return Object.keys(n).length > 0 && (t.extensions = n), t;
  }
}
class Eo extends ht {
  constructor() {
    super(), this.isCamera = !0, this.type = "Camera", this.matrixWorldInverse = new nt(), this.projectionMatrix = new nt(), this.projectionMatrixInverse = new nt(), this.coordinateSystem = qt;
  }
  copy(e, t) {
    return super.copy(e, t), this.matrixWorldInverse.copy(e.matrixWorldInverse), this.projectionMatrix.copy(e.projectionMatrix), this.projectionMatrixInverse.copy(e.projectionMatrixInverse), this.coordinateSystem = e.coordinateSystem, this;
  }
  getWorldDirection(e) {
    this.updateWorldMatrix(!0, !1);
    const t = this.matrixWorld.elements;
    return e.set(-t[8], -t[9], -t[10]).normalize();
  }
  updateMatrixWorld(e) {
    super.updateMatrixWorld(e), this.matrixWorldInverse.copy(this.matrixWorld).invert();
  }
  updateWorldMatrix(e, t) {
    super.updateWorldMatrix(e, t), this.matrixWorldInverse.copy(this.matrixWorld).invert();
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
class At extends Eo {
  constructor(e = 50, t = 1, n = 0.1, r = 2e3) {
    super(), this.isPerspectiveCamera = !0, this.type = "PerspectiveCamera", this.fov = e, this.zoom = 1, this.near = n, this.far = r, this.focus = 10, this.aspect = t, this.view = null, this.filmGauge = 35, this.filmOffset = 0, this.updateProjectionMatrix();
  }
  copy(e, t) {
    return super.copy(e, t), this.fov = e.fov, this.zoom = e.zoom, this.near = e.near, this.far = e.far, this.focus = e.focus, this.aspect = e.aspect, this.view = e.view === null ? null : Object.assign({}, e.view), this.filmGauge = e.filmGauge, this.filmOffset = e.filmOffset, this;
  }
  /**
   * Sets the FOV by focal length in respect to the current .filmGauge.
   *
   * The default film gauge is 35, so that the focal length can be specified for
   * a 35mm (full frame) camera.
   *
   * Values for focal length and film gauge must have the same unit.
   */
  setFocalLength(e) {
    const t = 0.5 * this.getFilmHeight() / e;
    this.fov = Mi * 2 * Math.atan(t), this.updateProjectionMatrix();
  }
  /**
   * Calculates the focal length from the current .fov and .filmGauge.
   */
  getFocalLength() {
    const e = Math.tan(ui * 0.5 * this.fov);
    return 0.5 * this.getFilmHeight() / e;
  }
  getEffectiveFOV() {
    return Mi * 2 * Math.atan(
      Math.tan(ui * 0.5 * this.fov) / this.zoom
    );
  }
  getFilmWidth() {
    return this.filmGauge * Math.min(this.aspect, 1);
  }
  getFilmHeight() {
    return this.filmGauge / Math.max(this.aspect, 1);
  }
  /**
   * Sets an offset in a larger frustum. This is useful for multi-window or
   * multi-monitor/multi-machine setups.
   *
   * For example, if you have 3x2 monitors and each monitor is 1920x1080 and
   * the monitors are in grid like this
   *
   *   +---+---+---+
   *   | A | B | C |
   *   +---+---+---+
   *   | D | E | F |
   *   +---+---+---+
   *
   * then for each monitor you would call it like this
   *
   *   const w = 1920;
   *   const h = 1080;
   *   const fullWidth = w * 3;
   *   const fullHeight = h * 2;
   *
   *   --A--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 0, h * 0, w, h );
   *   --B--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 1, h * 0, w, h );
   *   --C--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 2, h * 0, w, h );
   *   --D--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 0, h * 1, w, h );
   *   --E--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 1, h * 1, w, h );
   *   --F--
   *   camera.setViewOffset( fullWidth, fullHeight, w * 2, h * 1, w, h );
   *
   *   Note there is no reason monitors have to be the same size or in a grid.
   */
  setViewOffset(e, t, n, r, s, o) {
    this.aspect = e / t, this.view === null && (this.view = {
      enabled: !0,
      fullWidth: 1,
      fullHeight: 1,
      offsetX: 0,
      offsetY: 0,
      width: 1,
      height: 1
    }), this.view.enabled = !0, this.view.fullWidth = e, this.view.fullHeight = t, this.view.offsetX = n, this.view.offsetY = r, this.view.width = s, this.view.height = o, this.updateProjectionMatrix();
  }
  clearViewOffset() {
    this.view !== null && (this.view.enabled = !1), this.updateProjectionMatrix();
  }
  updateProjectionMatrix() {
    const e = this.near;
    let t = e * Math.tan(ui * 0.5 * this.fov) / this.zoom, n = 2 * t, r = this.aspect * n, s = -0.5 * r;
    const o = this.view;
    if (this.view !== null && this.view.enabled) {
      const l = o.fullWidth, c = o.fullHeight;
      s += o.offsetX * r / l, t -= o.offsetY * n / c, r *= o.width / l, n *= o.height / c;
    }
    const a = this.filmOffset;
    a !== 0 && (s += e * a / this.getFilmWidth()), this.projectionMatrix.makePerspective(s, s + r, t, t - n, e, this.far, this.coordinateSystem), this.projectionMatrixInverse.copy(this.projectionMatrix).invert();
  }
  toJSON(e) {
    const t = super.toJSON(e);
    return t.object.fov = this.fov, t.object.zoom = this.zoom, t.object.near = this.near, t.object.far = this.far, t.object.focus = this.focus, t.object.aspect = this.aspect, this.view !== null && (t.object.view = Object.assign({}, this.view)), t.object.filmGauge = this.filmGauge, t.object.filmOffset = this.filmOffset, t;
  }
}
const Xn = -90, Yn = 1;
class Lc extends ht {
  constructor(e, t, n) {
    super(), this.type = "CubeCamera", this.renderTarget = n, this.coordinateSystem = null;
    const r = new At(Xn, Yn, e, t);
    r.layers = this.layers, this.add(r);
    const s = new At(Xn, Yn, e, t);
    s.layers = this.layers, this.add(s);
    const o = new At(Xn, Yn, e, t);
    o.layers = this.layers, this.add(o);
    const a = new At(Xn, Yn, e, t);
    a.layers = this.layers, this.add(a);
    const l = new At(Xn, Yn, e, t);
    l.layers = this.layers, this.add(l);
    const c = new At(Xn, Yn, e, t);
    c.layers = this.layers, this.add(c);
  }
  updateCoordinateSystem() {
    const e = this.coordinateSystem, t = this.children.concat(), [n, r, s, o, a, l] = t;
    for (const c of t)
      this.remove(c);
    if (e === qt)
      n.up.set(0, 1, 0), n.lookAt(1, 0, 0), r.up.set(0, 1, 0), r.lookAt(-1, 0, 0), s.up.set(0, 0, -1), s.lookAt(0, 1, 0), o.up.set(0, 0, 1), o.lookAt(0, -1, 0), a.up.set(0, 1, 0), a.lookAt(0, 0, 1), l.up.set(0, 1, 0), l.lookAt(0, 0, -1);
    else if (e === tr)
      n.up.set(0, -1, 0), n.lookAt(-1, 0, 0), r.up.set(0, -1, 0), r.lookAt(1, 0, 0), s.up.set(0, 0, 1), s.lookAt(0, 1, 0), o.up.set(0, 0, -1), o.lookAt(0, -1, 0), a.up.set(0, -1, 0), a.lookAt(0, 0, 1), l.up.set(0, -1, 0), l.lookAt(0, 0, -1);
    else
      throw new Error("THREE.CubeCamera.updateCoordinateSystem(): Invalid coordinate system: " + e);
    for (const c of t)
      this.add(c), c.updateMatrixWorld();
  }
  update(e, t) {
    this.parent === null && this.updateMatrixWorld();
    const n = this.renderTarget;
    this.coordinateSystem !== e.coordinateSystem && (this.coordinateSystem = e.coordinateSystem, this.updateCoordinateSystem());
    const [r, s, o, a, l, c] = this.children, h = e.getRenderTarget(), f = e.xr.enabled;
    e.xr.enabled = !1;
    const u = n.texture.generateMipmaps;
    n.texture.generateMipmaps = !1, e.setRenderTarget(n, 0), e.render(t, r), e.setRenderTarget(n, 1), e.render(t, s), e.setRenderTarget(n, 2), e.render(t, o), e.setRenderTarget(n, 3), e.render(t, a), e.setRenderTarget(n, 4), e.render(t, l), n.texture.generateMipmaps = u, e.setRenderTarget(n, 5), e.render(t, c), e.setRenderTarget(h), e.xr.enabled = f, n.texture.needsPMREMUpdate = !0;
  }
}
class yo extends St {
  constructor(e, t, n, r, s, o, a, l, c, h) {
    e = e !== void 0 ? e : [], t = t !== void 0 ? t : ei, super(e, t, n, r, s, o, a, l, c, h), this.isCubeTexture = !0, this.flipY = !1;
  }
  get images() {
    return this.image;
  }
  set images(e) {
    this.image = e;
  }
}
class Uc extends Tn {
  constructor(e = 1, t = {}) {
    super(e, e, t), this.isWebGLCubeRenderTarget = !0;
    const n = { width: e, height: e, depth: 1 }, r = [n, n, n, n, n, n];
    t.encoding !== void 0 && (di("THREE.WebGLCubeRenderTarget: option.encoding has been replaced by option.colorSpace."), t.colorSpace = t.encoding === Sn ? Oe : En), this.texture = new yo(r, t.mapping, t.wrapS, t.wrapT, t.magFilter, t.minFilter, t.format, t.type, t.anisotropy, t.colorSpace), this.texture.isRenderTargetTexture = !0, this.texture.generateMipmaps = t.generateMipmaps !== void 0 ? t.generateMipmaps : !1, this.texture.minFilter = t.minFilter !== void 0 ? t.minFilter : Tt;
  }
  fromEquirectangularTexture(e, t) {
    this.texture.type = t.type, this.texture.colorSpace = t.colorSpace, this.texture.generateMipmaps = t.generateMipmaps, this.texture.minFilter = t.minFilter, this.texture.magFilter = t.magFilter;
    const n = {
      uniforms: {
        tEquirect: { value: null }
      },
      vertexShader: (
        /* glsl */
        `

				varying vec3 vWorldDirection;

				vec3 transformDirection( in vec3 dir, in mat4 matrix ) {

					return normalize( ( matrix * vec4( dir, 0.0 ) ).xyz );

				}

				void main() {

					vWorldDirection = transformDirection( position, modelMatrix );

					#include <begin_vertex>
					#include <project_vertex>

				}
			`
      ),
      fragmentShader: (
        /* glsl */
        `

				uniform sampler2D tEquirect;

				varying vec3 vWorldDirection;

				#include <common>

				void main() {

					vec3 direction = normalize( vWorldDirection );

					vec2 sampleUV = equirectUv( direction );

					gl_FragColor = texture2D( tEquirect, sampleUV );

				}
			`
      )
    }, r = new bi(5, 5, 5), s = new bn({
      name: "CubemapFromEquirect",
      uniforms: ii(n.uniforms),
      vertexShader: n.vertexShader,
      fragmentShader: n.fragmentShader,
      side: gt,
      blending: sn
    });
    s.uniforms.tEquirect.value = t;
    const o = new Nt(r, s), a = t.minFilter;
    return t.minFilter === vi && (t.minFilter = Tt), new Lc(1, 10, this).update(e, o), t.minFilter = a, o.geometry.dispose(), o.material.dispose(), this;
  }
  clear(e, t, n, r) {
    const s = e.getRenderTarget();
    for (let o = 0; o < 6; o++)
      e.setRenderTarget(this, o), e.clear(t, n, r);
    e.setRenderTarget(s);
  }
}
const Or = /* @__PURE__ */ new U(), Dc = /* @__PURE__ */ new U(), Ic = /* @__PURE__ */ new Be();
class en {
  constructor(e = new U(1, 0, 0), t = 0) {
    this.isPlane = !0, this.normal = e, this.constant = t;
  }
  set(e, t) {
    return this.normal.copy(e), this.constant = t, this;
  }
  setComponents(e, t, n, r) {
    return this.normal.set(e, t, n), this.constant = r, this;
  }
  setFromNormalAndCoplanarPoint(e, t) {
    return this.normal.copy(e), this.constant = -t.dot(this.normal), this;
  }
  setFromCoplanarPoints(e, t, n) {
    const r = Or.subVectors(n, t).cross(Dc.subVectors(e, t)).normalize();
    return this.setFromNormalAndCoplanarPoint(r, e), this;
  }
  copy(e) {
    return this.normal.copy(e.normal), this.constant = e.constant, this;
  }
  normalize() {
    const e = 1 / this.normal.length();
    return this.normal.multiplyScalar(e), this.constant *= e, this;
  }
  negate() {
    return this.constant *= -1, this.normal.negate(), this;
  }
  distanceToPoint(e) {
    return this.normal.dot(e) + this.constant;
  }
  distanceToSphere(e) {
    return this.distanceToPoint(e.center) - e.radius;
  }
  projectPoint(e, t) {
    return t.copy(e).addScaledVector(this.normal, -this.distanceToPoint(e));
  }
  intersectLine(e, t) {
    const n = e.delta(Or), r = this.normal.dot(n);
    if (r === 0)
      return this.distanceToPoint(e.start) === 0 ? t.copy(e.start) : null;
    const s = -(e.start.dot(this.normal) + this.constant) / r;
    return s < 0 || s > 1 ? null : t.copy(e.start).addScaledVector(n, s);
  }
  intersectsLine(e) {
    const t = this.distanceToPoint(e.start), n = this.distanceToPoint(e.end);
    return t < 0 && n > 0 || n < 0 && t > 0;
  }
  intersectsBox(e) {
    return e.intersectsPlane(this);
  }
  intersectsSphere(e) {
    return e.intersectsPlane(this);
  }
  coplanarPoint(e) {
    return e.copy(this.normal).multiplyScalar(-this.constant);
  }
  applyMatrix4(e, t) {
    const n = t || Ic.getNormalMatrix(e), r = this.coplanarPoint(Or).applyMatrix4(e), s = this.normal.applyMatrix3(n).normalize();
    return this.constant = -r.dot(s), this;
  }
  translate(e) {
    return this.constant -= e.dot(this.normal), this;
  }
  equals(e) {
    return e.normal.equals(this.normal) && e.constant === this.constant;
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
const mn = /* @__PURE__ */ new ls(), Yi = /* @__PURE__ */ new U();
class us {
  constructor(e = new en(), t = new en(), n = new en(), r = new en(), s = new en(), o = new en()) {
    this.planes = [e, t, n, r, s, o];
  }
  set(e, t, n, r, s, o) {
    const a = this.planes;
    return a[0].copy(e), a[1].copy(t), a[2].copy(n), a[3].copy(r), a[4].copy(s), a[5].copy(o), this;
  }
  copy(e) {
    const t = this.planes;
    for (let n = 0; n < 6; n++)
      t[n].copy(e.planes[n]);
    return this;
  }
  setFromProjectionMatrix(e, t = qt) {
    const n = this.planes, r = e.elements, s = r[0], o = r[1], a = r[2], l = r[3], c = r[4], h = r[5], f = r[6], u = r[7], m = r[8], g = r[9], x = r[10], p = r[11], d = r[12], A = r[13], _ = r[14], T = r[15];
    if (n[0].setComponents(l - s, u - c, p - m, T - d).normalize(), n[1].setComponents(l + s, u + c, p + m, T + d).normalize(), n[2].setComponents(l + o, u + h, p + g, T + A).normalize(), n[3].setComponents(l - o, u - h, p - g, T - A).normalize(), n[4].setComponents(l - a, u - f, p - x, T - _).normalize(), t === qt)
      n[5].setComponents(l + a, u + f, p + x, T + _).normalize();
    else if (t === tr)
      n[5].setComponents(a, f, x, _).normalize();
    else
      throw new Error("THREE.Frustum.setFromProjectionMatrix(): Invalid coordinate system: " + t);
    return this;
  }
  intersectsObject(e) {
    if (e.boundingSphere !== void 0)
      e.boundingSphere === null && e.computeBoundingSphere(), mn.copy(e.boundingSphere).applyMatrix4(e.matrixWorld);
    else {
      const t = e.geometry;
      t.boundingSphere === null && t.computeBoundingSphere(), mn.copy(t.boundingSphere).applyMatrix4(e.matrixWorld);
    }
    return this.intersectsSphere(mn);
  }
  intersectsSprite(e) {
    return mn.center.set(0, 0, 0), mn.radius = 0.7071067811865476, mn.applyMatrix4(e.matrixWorld), this.intersectsSphere(mn);
  }
  intersectsSphere(e) {
    const t = this.planes, n = e.center, r = -e.radius;
    for (let s = 0; s < 6; s++)
      if (t[s].distanceToPoint(n) < r)
        return !1;
    return !0;
  }
  intersectsBox(e) {
    const t = this.planes;
    for (let n = 0; n < 6; n++) {
      const r = t[n];
      if (Yi.x = r.normal.x > 0 ? e.max.x : e.min.x, Yi.y = r.normal.y > 0 ? e.max.y : e.min.y, Yi.z = r.normal.z > 0 ? e.max.z : e.min.z, r.distanceToPoint(Yi) < 0)
        return !1;
    }
    return !0;
  }
  containsPoint(e) {
    const t = this.planes;
    for (let n = 0; n < 6; n++)
      if (t[n].distanceToPoint(e) < 0)
        return !1;
    return !0;
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
function To() {
  let i = null, e = !1, t = null, n = null;
  function r(s, o) {
    t(s, o), n = i.requestAnimationFrame(r);
  }
  return {
    start: function() {
      e !== !0 && t !== null && (n = i.requestAnimationFrame(r), e = !0);
    },
    stop: function() {
      i.cancelAnimationFrame(n), e = !1;
    },
    setAnimationLoop: function(s) {
      t = s;
    },
    setContext: function(s) {
      i = s;
    }
  };
}
function Nc(i, e) {
  const t = e.isWebGL2, n = /* @__PURE__ */ new WeakMap();
  function r(c, h) {
    const f = c.array, u = c.usage, m = i.createBuffer();
    i.bindBuffer(h, m), i.bufferData(h, f, u), c.onUploadCallback();
    let g;
    if (f instanceof Float32Array)
      g = i.FLOAT;
    else if (f instanceof Uint16Array)
      if (c.isFloat16BufferAttribute)
        if (t)
          g = i.HALF_FLOAT;
        else
          throw new Error("THREE.WebGLAttributes: Usage of Float16BufferAttribute requires WebGL2.");
      else
        g = i.UNSIGNED_SHORT;
    else if (f instanceof Int16Array)
      g = i.SHORT;
    else if (f instanceof Uint32Array)
      g = i.UNSIGNED_INT;
    else if (f instanceof Int32Array)
      g = i.INT;
    else if (f instanceof Int8Array)
      g = i.BYTE;
    else if (f instanceof Uint8Array)
      g = i.UNSIGNED_BYTE;
    else if (f instanceof Uint8ClampedArray)
      g = i.UNSIGNED_BYTE;
    else
      throw new Error("THREE.WebGLAttributes: Unsupported buffer data format: " + f);
    return {
      buffer: m,
      type: g,
      bytesPerElement: f.BYTES_PER_ELEMENT,
      version: c.version
    };
  }
  function s(c, h, f) {
    const u = h.array, m = h.updateRange;
    i.bindBuffer(f, c), m.count === -1 ? i.bufferSubData(f, 0, u) : (t ? i.bufferSubData(
      f,
      m.offset * u.BYTES_PER_ELEMENT,
      u,
      m.offset,
      m.count
    ) : i.bufferSubData(
      f,
      m.offset * u.BYTES_PER_ELEMENT,
      u.subarray(m.offset, m.offset + m.count)
    ), m.count = -1), h.onUploadCallback();
  }
  function o(c) {
    return c.isInterleavedBufferAttribute && (c = c.data), n.get(c);
  }
  function a(c) {
    c.isInterleavedBufferAttribute && (c = c.data);
    const h = n.get(c);
    h && (i.deleteBuffer(h.buffer), n.delete(c));
  }
  function l(c, h) {
    if (c.isGLBufferAttribute) {
      const u = n.get(c);
      (!u || u.version < c.version) && n.set(c, {
        buffer: c.buffer,
        type: c.type,
        bytesPerElement: c.elementSize,
        version: c.version
      });
      return;
    }
    c.isInterleavedBufferAttribute && (c = c.data);
    const f = n.get(c);
    f === void 0 ? n.set(c, r(c, h)) : f.version < c.version && (s(f.buffer, c, h), f.version = c.version);
  }
  return {
    get: o,
    remove: a,
    update: l
  };
}
class ar extends cn {
  constructor(e = 1, t = 1, n = 1, r = 1) {
    super(), this.type = "PlaneGeometry", this.parameters = {
      width: e,
      height: t,
      widthSegments: n,
      heightSegments: r
    };
    const s = e / 2, o = t / 2, a = Math.floor(n), l = Math.floor(r), c = a + 1, h = l + 1, f = e / a, u = t / l, m = [], g = [], x = [], p = [];
    for (let d = 0; d < h; d++) {
      const A = d * u - o;
      for (let _ = 0; _ < c; _++) {
        const T = _ * f - s;
        g.push(T, -A, 0), x.push(0, 0, 1), p.push(_ / a), p.push(1 - d / l);
      }
    }
    for (let d = 0; d < l; d++)
      for (let A = 0; A < a; A++) {
        const _ = A + c * d, T = A + c * (d + 1), C = A + 1 + c * (d + 1), L = A + 1 + c * d;
        m.push(_, T, L), m.push(T, C, L);
      }
    this.setIndex(m), this.setAttribute("position", new jt(g, 3)), this.setAttribute("normal", new jt(x, 3)), this.setAttribute("uv", new jt(p, 2));
  }
  copy(e) {
    return super.copy(e), this.parameters = Object.assign({}, e.parameters), this;
  }
  static fromJSON(e) {
    return new ar(e.width, e.height, e.widthSegments, e.heightSegments);
  }
}
var Oc = `#ifdef USE_ALPHAHASH
	if ( diffuseColor.a < getAlphaHashThreshold( vPosition ) ) discard;
#endif`, Fc = `#ifdef USE_ALPHAHASH
	const float ALPHA_HASH_SCALE = 0.05;
	float hash2D( vec2 value ) {
		return fract( 1.0e4 * sin( 17.0 * value.x + 0.1 * value.y ) * ( 0.1 + abs( sin( 13.0 * value.y + value.x ) ) ) );
	}
	float hash3D( vec3 value ) {
		return hash2D( vec2( hash2D( value.xy ), value.z ) );
	}
	float getAlphaHashThreshold( vec3 position ) {
		float maxDeriv = max(
			length( dFdx( position.xyz ) ),
			length( dFdy( position.xyz ) )
		);
		float pixScale = 1.0 / ( ALPHA_HASH_SCALE * maxDeriv );
		vec2 pixScales = vec2(
			exp2( floor( log2( pixScale ) ) ),
			exp2( ceil( log2( pixScale ) ) )
		);
		vec2 alpha = vec2(
			hash3D( floor( pixScales.x * position.xyz ) ),
			hash3D( floor( pixScales.y * position.xyz ) )
		);
		float lerpFactor = fract( log2( pixScale ) );
		float x = ( 1.0 - lerpFactor ) * alpha.x + lerpFactor * alpha.y;
		float a = min( lerpFactor, 1.0 - lerpFactor );
		vec3 cases = vec3(
			x * x / ( 2.0 * a * ( 1.0 - a ) ),
			( x - 0.5 * a ) / ( 1.0 - a ),
			1.0 - ( ( 1.0 - x ) * ( 1.0 - x ) / ( 2.0 * a * ( 1.0 - a ) ) )
		);
		float threshold = ( x < ( 1.0 - a ) )
			? ( ( x < a ) ? cases.x : cases.y )
			: cases.z;
		return clamp( threshold , 1.0e-6, 1.0 );
	}
#endif`, Bc = `#ifdef USE_ALPHAMAP
	diffuseColor.a *= texture2D( alphaMap, vAlphaMapUv ).g;
#endif`, zc = `#ifdef USE_ALPHAMAP
	uniform sampler2D alphaMap;
#endif`, Hc = `#ifdef USE_ALPHATEST
	if ( diffuseColor.a < alphaTest ) discard;
#endif`, Gc = `#ifdef USE_ALPHATEST
	uniform float alphaTest;
#endif`, Vc = `#ifdef USE_AOMAP
	float ambientOcclusion = ( texture2D( aoMap, vAoMapUv ).r - 1.0 ) * aoMapIntensity + 1.0;
	reflectedLight.indirectDiffuse *= ambientOcclusion;
	#if defined( USE_ENVMAP ) && defined( STANDARD )
		float dotNV = saturate( dot( geometry.normal, geometry.viewDir ) );
		reflectedLight.indirectSpecular *= computeSpecularOcclusion( dotNV, ambientOcclusion, material.roughness );
	#endif
#endif`, kc = `#ifdef USE_AOMAP
	uniform sampler2D aoMap;
	uniform float aoMapIntensity;
#endif`, Wc = `vec3 transformed = vec3( position );
#ifdef USE_ALPHAHASH
	vPosition = vec3( position );
#endif`, Xc = `vec3 objectNormal = vec3( normal );
#ifdef USE_TANGENT
	vec3 objectTangent = vec3( tangent.xyz );
#endif`, Yc = `float G_BlinnPhong_Implicit( ) {
	return 0.25;
}
float D_BlinnPhong( const in float shininess, const in float dotNH ) {
	return RECIPROCAL_PI * ( shininess * 0.5 + 1.0 ) * pow( dotNH, shininess );
}
vec3 BRDF_BlinnPhong( const in vec3 lightDir, const in vec3 viewDir, const in vec3 normal, const in vec3 specularColor, const in float shininess ) {
	vec3 halfDir = normalize( lightDir + viewDir );
	float dotNH = saturate( dot( normal, halfDir ) );
	float dotVH = saturate( dot( viewDir, halfDir ) );
	vec3 F = F_Schlick( specularColor, 1.0, dotVH );
	float G = G_BlinnPhong_Implicit( );
	float D = D_BlinnPhong( shininess, dotNH );
	return F * ( G * D );
} // validated`, qc = `#ifdef USE_IRIDESCENCE
	const mat3 XYZ_TO_REC709 = mat3(
		 3.2404542, -0.9692660,  0.0556434,
		-1.5371385,  1.8760108, -0.2040259,
		-0.4985314,  0.0415560,  1.0572252
	);
	vec3 Fresnel0ToIor( vec3 fresnel0 ) {
		vec3 sqrtF0 = sqrt( fresnel0 );
		return ( vec3( 1.0 ) + sqrtF0 ) / ( vec3( 1.0 ) - sqrtF0 );
	}
	vec3 IorToFresnel0( vec3 transmittedIor, float incidentIor ) {
		return pow2( ( transmittedIor - vec3( incidentIor ) ) / ( transmittedIor + vec3( incidentIor ) ) );
	}
	float IorToFresnel0( float transmittedIor, float incidentIor ) {
		return pow2( ( transmittedIor - incidentIor ) / ( transmittedIor + incidentIor ));
	}
	vec3 evalSensitivity( float OPD, vec3 shift ) {
		float phase = 2.0 * PI * OPD * 1.0e-9;
		vec3 val = vec3( 5.4856e-13, 4.4201e-13, 5.2481e-13 );
		vec3 pos = vec3( 1.6810e+06, 1.7953e+06, 2.2084e+06 );
		vec3 var = vec3( 4.3278e+09, 9.3046e+09, 6.6121e+09 );
		vec3 xyz = val * sqrt( 2.0 * PI * var ) * cos( pos * phase + shift ) * exp( - pow2( phase ) * var );
		xyz.x += 9.7470e-14 * sqrt( 2.0 * PI * 4.5282e+09 ) * cos( 2.2399e+06 * phase + shift[ 0 ] ) * exp( - 4.5282e+09 * pow2( phase ) );
		xyz /= 1.0685e-7;
		vec3 rgb = XYZ_TO_REC709 * xyz;
		return rgb;
	}
	vec3 evalIridescence( float outsideIOR, float eta2, float cosTheta1, float thinFilmThickness, vec3 baseF0 ) {
		vec3 I;
		float iridescenceIOR = mix( outsideIOR, eta2, smoothstep( 0.0, 0.03, thinFilmThickness ) );
		float sinTheta2Sq = pow2( outsideIOR / iridescenceIOR ) * ( 1.0 - pow2( cosTheta1 ) );
		float cosTheta2Sq = 1.0 - sinTheta2Sq;
		if ( cosTheta2Sq < 0.0 ) {
			return vec3( 1.0 );
		}
		float cosTheta2 = sqrt( cosTheta2Sq );
		float R0 = IorToFresnel0( iridescenceIOR, outsideIOR );
		float R12 = F_Schlick( R0, 1.0, cosTheta1 );
		float T121 = 1.0 - R12;
		float phi12 = 0.0;
		if ( iridescenceIOR < outsideIOR ) phi12 = PI;
		float phi21 = PI - phi12;
		vec3 baseIOR = Fresnel0ToIor( clamp( baseF0, 0.0, 0.9999 ) );		vec3 R1 = IorToFresnel0( baseIOR, iridescenceIOR );
		vec3 R23 = F_Schlick( R1, 1.0, cosTheta2 );
		vec3 phi23 = vec3( 0.0 );
		if ( baseIOR[ 0 ] < iridescenceIOR ) phi23[ 0 ] = PI;
		if ( baseIOR[ 1 ] < iridescenceIOR ) phi23[ 1 ] = PI;
		if ( baseIOR[ 2 ] < iridescenceIOR ) phi23[ 2 ] = PI;
		float OPD = 2.0 * iridescenceIOR * thinFilmThickness * cosTheta2;
		vec3 phi = vec3( phi21 ) + phi23;
		vec3 R123 = clamp( R12 * R23, 1e-5, 0.9999 );
		vec3 r123 = sqrt( R123 );
		vec3 Rs = pow2( T121 ) * R23 / ( vec3( 1.0 ) - R123 );
		vec3 C0 = R12 + Rs;
		I = C0;
		vec3 Cm = Rs - T121;
		for ( int m = 1; m <= 2; ++ m ) {
			Cm *= r123;
			vec3 Sm = 2.0 * evalSensitivity( float( m ) * OPD, float( m ) * phi );
			I += Cm * Sm;
		}
		return max( I, vec3( 0.0 ) );
	}
#endif`, jc = `#ifdef USE_BUMPMAP
	uniform sampler2D bumpMap;
	uniform float bumpScale;
	vec2 dHdxy_fwd() {
		vec2 dSTdx = dFdx( vBumpMapUv );
		vec2 dSTdy = dFdy( vBumpMapUv );
		float Hll = bumpScale * texture2D( bumpMap, vBumpMapUv ).x;
		float dBx = bumpScale * texture2D( bumpMap, vBumpMapUv + dSTdx ).x - Hll;
		float dBy = bumpScale * texture2D( bumpMap, vBumpMapUv + dSTdy ).x - Hll;
		return vec2( dBx, dBy );
	}
	vec3 perturbNormalArb( vec3 surf_pos, vec3 surf_norm, vec2 dHdxy, float faceDirection ) {
		vec3 vSigmaX = dFdx( surf_pos.xyz );
		vec3 vSigmaY = dFdy( surf_pos.xyz );
		vec3 vN = surf_norm;
		vec3 R1 = cross( vSigmaY, vN );
		vec3 R2 = cross( vN, vSigmaX );
		float fDet = dot( vSigmaX, R1 ) * faceDirection;
		vec3 vGrad = sign( fDet ) * ( dHdxy.x * R1 + dHdxy.y * R2 );
		return normalize( abs( fDet ) * surf_norm - vGrad );
	}
#endif`, Zc = `#if NUM_CLIPPING_PLANES > 0
	vec4 plane;
	#pragma unroll_loop_start
	for ( int i = 0; i < UNION_CLIPPING_PLANES; i ++ ) {
		plane = clippingPlanes[ i ];
		if ( dot( vClipPosition, plane.xyz ) > plane.w ) discard;
	}
	#pragma unroll_loop_end
	#if UNION_CLIPPING_PLANES < NUM_CLIPPING_PLANES
		bool clipped = true;
		#pragma unroll_loop_start
		for ( int i = UNION_CLIPPING_PLANES; i < NUM_CLIPPING_PLANES; i ++ ) {
			plane = clippingPlanes[ i ];
			clipped = ( dot( vClipPosition, plane.xyz ) > plane.w ) && clipped;
		}
		#pragma unroll_loop_end
		if ( clipped ) discard;
	#endif
#endif`, Kc = `#if NUM_CLIPPING_PLANES > 0
	varying vec3 vClipPosition;
	uniform vec4 clippingPlanes[ NUM_CLIPPING_PLANES ];
#endif`, Jc = `#if NUM_CLIPPING_PLANES > 0
	varying vec3 vClipPosition;
#endif`, $c = `#if NUM_CLIPPING_PLANES > 0
	vClipPosition = - mvPosition.xyz;
#endif`, Qc = `#if defined( USE_COLOR_ALPHA )
	diffuseColor *= vColor;
#elif defined( USE_COLOR )
	diffuseColor.rgb *= vColor;
#endif`, eh = `#if defined( USE_COLOR_ALPHA )
	varying vec4 vColor;
#elif defined( USE_COLOR )
	varying vec3 vColor;
#endif`, th = `#if defined( USE_COLOR_ALPHA )
	varying vec4 vColor;
#elif defined( USE_COLOR ) || defined( USE_INSTANCING_COLOR )
	varying vec3 vColor;
#endif`, nh = `#if defined( USE_COLOR_ALPHA )
	vColor = vec4( 1.0 );
#elif defined( USE_COLOR ) || defined( USE_INSTANCING_COLOR )
	vColor = vec3( 1.0 );
#endif
#ifdef USE_COLOR
	vColor *= color;
#endif
#ifdef USE_INSTANCING_COLOR
	vColor.xyz *= instanceColor.xyz;
#endif`, ih = `#define PI 3.141592653589793
#define PI2 6.283185307179586
#define PI_HALF 1.5707963267948966
#define RECIPROCAL_PI 0.3183098861837907
#define RECIPROCAL_PI2 0.15915494309189535
#define EPSILON 1e-6
#ifndef saturate
#define saturate( a ) clamp( a, 0.0, 1.0 )
#endif
#define whiteComplement( a ) ( 1.0 - saturate( a ) )
float pow2( const in float x ) { return x*x; }
vec3 pow2( const in vec3 x ) { return x*x; }
float pow3( const in float x ) { return x*x*x; }
float pow4( const in float x ) { float x2 = x*x; return x2*x2; }
float max3( const in vec3 v ) { return max( max( v.x, v.y ), v.z ); }
float average( const in vec3 v ) { return dot( v, vec3( 0.3333333 ) ); }
highp float rand( const in vec2 uv ) {
	const highp float a = 12.9898, b = 78.233, c = 43758.5453;
	highp float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );
	return fract( sin( sn ) * c );
}
#ifdef HIGH_PRECISION
	float precisionSafeLength( vec3 v ) { return length( v ); }
#else
	float precisionSafeLength( vec3 v ) {
		float maxComponent = max3( abs( v ) );
		return length( v / maxComponent ) * maxComponent;
	}
#endif
struct IncidentLight {
	vec3 color;
	vec3 direction;
	bool visible;
};
struct ReflectedLight {
	vec3 directDiffuse;
	vec3 directSpecular;
	vec3 indirectDiffuse;
	vec3 indirectSpecular;
};
struct GeometricContext {
	vec3 position;
	vec3 normal;
	vec3 viewDir;
#ifdef USE_CLEARCOAT
	vec3 clearcoatNormal;
#endif
};
#ifdef USE_ALPHAHASH
	varying vec3 vPosition;
#endif
vec3 transformDirection( in vec3 dir, in mat4 matrix ) {
	return normalize( ( matrix * vec4( dir, 0.0 ) ).xyz );
}
vec3 inverseTransformDirection( in vec3 dir, in mat4 matrix ) {
	return normalize( ( vec4( dir, 0.0 ) * matrix ).xyz );
}
mat3 transposeMat3( const in mat3 m ) {
	mat3 tmp;
	tmp[ 0 ] = vec3( m[ 0 ].x, m[ 1 ].x, m[ 2 ].x );
	tmp[ 1 ] = vec3( m[ 0 ].y, m[ 1 ].y, m[ 2 ].y );
	tmp[ 2 ] = vec3( m[ 0 ].z, m[ 1 ].z, m[ 2 ].z );
	return tmp;
}
float luminance( const in vec3 rgb ) {
	const vec3 weights = vec3( 0.2126729, 0.7151522, 0.0721750 );
	return dot( weights, rgb );
}
bool isPerspectiveMatrix( mat4 m ) {
	return m[ 2 ][ 3 ] == - 1.0;
}
vec2 equirectUv( in vec3 dir ) {
	float u = atan( dir.z, dir.x ) * RECIPROCAL_PI2 + 0.5;
	float v = asin( clamp( dir.y, - 1.0, 1.0 ) ) * RECIPROCAL_PI + 0.5;
	return vec2( u, v );
}
vec3 BRDF_Lambert( const in vec3 diffuseColor ) {
	return RECIPROCAL_PI * diffuseColor;
}
vec3 F_Schlick( const in vec3 f0, const in float f90, const in float dotVH ) {
	float fresnel = exp2( ( - 5.55473 * dotVH - 6.98316 ) * dotVH );
	return f0 * ( 1.0 - fresnel ) + ( f90 * fresnel );
}
float F_Schlick( const in float f0, const in float f90, const in float dotVH ) {
	float fresnel = exp2( ( - 5.55473 * dotVH - 6.98316 ) * dotVH );
	return f0 * ( 1.0 - fresnel ) + ( f90 * fresnel );
} // validated`, rh = `#ifdef ENVMAP_TYPE_CUBE_UV
	#define cubeUV_minMipLevel 4.0
	#define cubeUV_minTileSize 16.0
	float getFace( vec3 direction ) {
		vec3 absDirection = abs( direction );
		float face = - 1.0;
		if ( absDirection.x > absDirection.z ) {
			if ( absDirection.x > absDirection.y )
				face = direction.x > 0.0 ? 0.0 : 3.0;
			else
				face = direction.y > 0.0 ? 1.0 : 4.0;
		} else {
			if ( absDirection.z > absDirection.y )
				face = direction.z > 0.0 ? 2.0 : 5.0;
			else
				face = direction.y > 0.0 ? 1.0 : 4.0;
		}
		return face;
	}
	vec2 getUV( vec3 direction, float face ) {
		vec2 uv;
		if ( face == 0.0 ) {
			uv = vec2( direction.z, direction.y ) / abs( direction.x );
		} else if ( face == 1.0 ) {
			uv = vec2( - direction.x, - direction.z ) / abs( direction.y );
		} else if ( face == 2.0 ) {
			uv = vec2( - direction.x, direction.y ) / abs( direction.z );
		} else if ( face == 3.0 ) {
			uv = vec2( - direction.z, direction.y ) / abs( direction.x );
		} else if ( face == 4.0 ) {
			uv = vec2( - direction.x, direction.z ) / abs( direction.y );
		} else {
			uv = vec2( direction.x, direction.y ) / abs( direction.z );
		}
		return 0.5 * ( uv + 1.0 );
	}
	vec3 bilinearCubeUV( sampler2D envMap, vec3 direction, float mipInt ) {
		float face = getFace( direction );
		float filterInt = max( cubeUV_minMipLevel - mipInt, 0.0 );
		mipInt = max( mipInt, cubeUV_minMipLevel );
		float faceSize = exp2( mipInt );
		highp vec2 uv = getUV( direction, face ) * ( faceSize - 2.0 ) + 1.0;
		if ( face > 2.0 ) {
			uv.y += faceSize;
			face -= 3.0;
		}
		uv.x += face * faceSize;
		uv.x += filterInt * 3.0 * cubeUV_minTileSize;
		uv.y += 4.0 * ( exp2( CUBEUV_MAX_MIP ) - faceSize );
		uv.x *= CUBEUV_TEXEL_WIDTH;
		uv.y *= CUBEUV_TEXEL_HEIGHT;
		#ifdef texture2DGradEXT
			return texture2DGradEXT( envMap, uv, vec2( 0.0 ), vec2( 0.0 ) ).rgb;
		#else
			return texture2D( envMap, uv ).rgb;
		#endif
	}
	#define cubeUV_r0 1.0
	#define cubeUV_v0 0.339
	#define cubeUV_m0 - 2.0
	#define cubeUV_r1 0.8
	#define cubeUV_v1 0.276
	#define cubeUV_m1 - 1.0
	#define cubeUV_r4 0.4
	#define cubeUV_v4 0.046
	#define cubeUV_m4 2.0
	#define cubeUV_r5 0.305
	#define cubeUV_v5 0.016
	#define cubeUV_m5 3.0
	#define cubeUV_r6 0.21
	#define cubeUV_v6 0.0038
	#define cubeUV_m6 4.0
	float roughnessToMip( float roughness ) {
		float mip = 0.0;
		if ( roughness >= cubeUV_r1 ) {
			mip = ( cubeUV_r0 - roughness ) * ( cubeUV_m1 - cubeUV_m0 ) / ( cubeUV_r0 - cubeUV_r1 ) + cubeUV_m0;
		} else if ( roughness >= cubeUV_r4 ) {
			mip = ( cubeUV_r1 - roughness ) * ( cubeUV_m4 - cubeUV_m1 ) / ( cubeUV_r1 - cubeUV_r4 ) + cubeUV_m1;
		} else if ( roughness >= cubeUV_r5 ) {
			mip = ( cubeUV_r4 - roughness ) * ( cubeUV_m5 - cubeUV_m4 ) / ( cubeUV_r4 - cubeUV_r5 ) + cubeUV_m4;
		} else if ( roughness >= cubeUV_r6 ) {
			mip = ( cubeUV_r5 - roughness ) * ( cubeUV_m6 - cubeUV_m5 ) / ( cubeUV_r5 - cubeUV_r6 ) + cubeUV_m5;
		} else {
			mip = - 2.0 * log2( 1.16 * roughness );		}
		return mip;
	}
	vec4 textureCubeUV( sampler2D envMap, vec3 sampleDir, float roughness ) {
		float mip = clamp( roughnessToMip( roughness ), cubeUV_m0, CUBEUV_MAX_MIP );
		float mipF = fract( mip );
		float mipInt = floor( mip );
		vec3 color0 = bilinearCubeUV( envMap, sampleDir, mipInt );
		if ( mipF == 0.0 ) {
			return vec4( color0, 1.0 );
		} else {
			vec3 color1 = bilinearCubeUV( envMap, sampleDir, mipInt + 1.0 );
			return vec4( mix( color0, color1, mipF ), 1.0 );
		}
	}
#endif`, sh = `vec3 transformedNormal = objectNormal;
#ifdef USE_INSTANCING
	mat3 m = mat3( instanceMatrix );
	transformedNormal /= vec3( dot( m[ 0 ], m[ 0 ] ), dot( m[ 1 ], m[ 1 ] ), dot( m[ 2 ], m[ 2 ] ) );
	transformedNormal = m * transformedNormal;
#endif
transformedNormal = normalMatrix * transformedNormal;
#ifdef FLIP_SIDED
	transformedNormal = - transformedNormal;
#endif
#ifdef USE_TANGENT
	vec3 transformedTangent = ( modelViewMatrix * vec4( objectTangent, 0.0 ) ).xyz;
	#ifdef FLIP_SIDED
		transformedTangent = - transformedTangent;
	#endif
#endif`, ah = `#ifdef USE_DISPLACEMENTMAP
	uniform sampler2D displacementMap;
	uniform float displacementScale;
	uniform float displacementBias;
#endif`, oh = `#ifdef USE_DISPLACEMENTMAP
	transformed += normalize( objectNormal ) * ( texture2D( displacementMap, vDisplacementMapUv ).x * displacementScale + displacementBias );
#endif`, lh = `#ifdef USE_EMISSIVEMAP
	vec4 emissiveColor = texture2D( emissiveMap, vEmissiveMapUv );
	totalEmissiveRadiance *= emissiveColor.rgb;
#endif`, ch = `#ifdef USE_EMISSIVEMAP
	uniform sampler2D emissiveMap;
#endif`, hh = "gl_FragColor = linearToOutputTexel( gl_FragColor );", uh = `vec4 LinearToLinear( in vec4 value ) {
	return value;
}
vec4 LinearTosRGB( in vec4 value ) {
	return vec4( mix( pow( value.rgb, vec3( 0.41666 ) ) * 1.055 - vec3( 0.055 ), value.rgb * 12.92, vec3( lessThanEqual( value.rgb, vec3( 0.0031308 ) ) ) ), value.a );
}`, fh = `#ifdef USE_ENVMAP
	#ifdef ENV_WORLDPOS
		vec3 cameraToFrag;
		if ( isOrthographic ) {
			cameraToFrag = normalize( vec3( - viewMatrix[ 0 ][ 2 ], - viewMatrix[ 1 ][ 2 ], - viewMatrix[ 2 ][ 2 ] ) );
		} else {
			cameraToFrag = normalize( vWorldPosition - cameraPosition );
		}
		vec3 worldNormal = inverseTransformDirection( normal, viewMatrix );
		#ifdef ENVMAP_MODE_REFLECTION
			vec3 reflectVec = reflect( cameraToFrag, worldNormal );
		#else
			vec3 reflectVec = refract( cameraToFrag, worldNormal, refractionRatio );
		#endif
	#else
		vec3 reflectVec = vReflect;
	#endif
	#ifdef ENVMAP_TYPE_CUBE
		vec4 envColor = textureCube( envMap, vec3( flipEnvMap * reflectVec.x, reflectVec.yz ) );
	#else
		vec4 envColor = vec4( 0.0 );
	#endif
	#ifdef ENVMAP_BLENDING_MULTIPLY
		outgoingLight = mix( outgoingLight, outgoingLight * envColor.xyz, specularStrength * reflectivity );
	#elif defined( ENVMAP_BLENDING_MIX )
		outgoingLight = mix( outgoingLight, envColor.xyz, specularStrength * reflectivity );
	#elif defined( ENVMAP_BLENDING_ADD )
		outgoingLight += envColor.xyz * specularStrength * reflectivity;
	#endif
#endif`, dh = `#ifdef USE_ENVMAP
	uniform float envMapIntensity;
	uniform float flipEnvMap;
	#ifdef ENVMAP_TYPE_CUBE
		uniform samplerCube envMap;
	#else
		uniform sampler2D envMap;
	#endif
	
#endif`, ph = `#ifdef USE_ENVMAP
	uniform float reflectivity;
	#if defined( USE_BUMPMAP ) || defined( USE_NORMALMAP ) || defined( PHONG ) || defined( LAMBERT )
		#define ENV_WORLDPOS
	#endif
	#ifdef ENV_WORLDPOS
		varying vec3 vWorldPosition;
		uniform float refractionRatio;
	#else
		varying vec3 vReflect;
	#endif
#endif`, mh = `#ifdef USE_ENVMAP
	#if defined( USE_BUMPMAP ) || defined( USE_NORMALMAP ) || defined( PHONG ) || defined( LAMBERT )
		#define ENV_WORLDPOS
	#endif
	#ifdef ENV_WORLDPOS
		
		varying vec3 vWorldPosition;
	#else
		varying vec3 vReflect;
		uniform float refractionRatio;
	#endif
#endif`, gh = `#ifdef USE_ENVMAP
	#ifdef ENV_WORLDPOS
		vWorldPosition = worldPosition.xyz;
	#else
		vec3 cameraToVertex;
		if ( isOrthographic ) {
			cameraToVertex = normalize( vec3( - viewMatrix[ 0 ][ 2 ], - viewMatrix[ 1 ][ 2 ], - viewMatrix[ 2 ][ 2 ] ) );
		} else {
			cameraToVertex = normalize( worldPosition.xyz - cameraPosition );
		}
		vec3 worldNormal = inverseTransformDirection( transformedNormal, viewMatrix );
		#ifdef ENVMAP_MODE_REFLECTION
			vReflect = reflect( cameraToVertex, worldNormal );
		#else
			vReflect = refract( cameraToVertex, worldNormal, refractionRatio );
		#endif
	#endif
#endif`, _h = `#ifdef USE_FOG
	vFogDepth = - mvPosition.z;
#endif`, vh = `#ifdef USE_FOG
	varying float vFogDepth;
#endif`, xh = `#ifdef USE_FOG
	#ifdef FOG_EXP2
		float fogFactor = 1.0 - exp( - fogDensity * fogDensity * vFogDepth * vFogDepth );
	#else
		float fogFactor = smoothstep( fogNear, fogFar, vFogDepth );
	#endif
	gl_FragColor.rgb = mix( gl_FragColor.rgb, fogColor, fogFactor );
#endif`, Mh = `#ifdef USE_FOG
	uniform vec3 fogColor;
	varying float vFogDepth;
	#ifdef FOG_EXP2
		uniform float fogDensity;
	#else
		uniform float fogNear;
		uniform float fogFar;
	#endif
#endif`, Sh = `#ifdef USE_GRADIENTMAP
	uniform sampler2D gradientMap;
#endif
vec3 getGradientIrradiance( vec3 normal, vec3 lightDirection ) {
	float dotNL = dot( normal, lightDirection );
	vec2 coord = vec2( dotNL * 0.5 + 0.5, 0.0 );
	#ifdef USE_GRADIENTMAP
		return vec3( texture2D( gradientMap, coord ).r );
	#else
		vec2 fw = fwidth( coord ) * 0.5;
		return mix( vec3( 0.7 ), vec3( 1.0 ), smoothstep( 0.7 - fw.x, 0.7 + fw.x, coord.x ) );
	#endif
}`, Eh = `#ifdef USE_LIGHTMAP
	vec4 lightMapTexel = texture2D( lightMap, vLightMapUv );
	vec3 lightMapIrradiance = lightMapTexel.rgb * lightMapIntensity;
	reflectedLight.indirectDiffuse += lightMapIrradiance;
#endif`, yh = `#ifdef USE_LIGHTMAP
	uniform sampler2D lightMap;
	uniform float lightMapIntensity;
#endif`, Th = `LambertMaterial material;
material.diffuseColor = diffuseColor.rgb;
material.specularStrength = specularStrength;`, Ah = `varying vec3 vViewPosition;
struct LambertMaterial {
	vec3 diffuseColor;
	float specularStrength;
};
void RE_Direct_Lambert( const in IncidentLight directLight, const in GeometricContext geometry, const in LambertMaterial material, inout ReflectedLight reflectedLight ) {
	float dotNL = saturate( dot( geometry.normal, directLight.direction ) );
	vec3 irradiance = dotNL * directLight.color;
	reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
void RE_IndirectDiffuse_Lambert( const in vec3 irradiance, const in GeometricContext geometry, const in LambertMaterial material, inout ReflectedLight reflectedLight ) {
	reflectedLight.indirectDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
#define RE_Direct				RE_Direct_Lambert
#define RE_IndirectDiffuse		RE_IndirectDiffuse_Lambert`, bh = `uniform bool receiveShadow;
uniform vec3 ambientLightColor;
uniform vec3 lightProbe[ 9 ];
vec3 shGetIrradianceAt( in vec3 normal, in vec3 shCoefficients[ 9 ] ) {
	float x = normal.x, y = normal.y, z = normal.z;
	vec3 result = shCoefficients[ 0 ] * 0.886227;
	result += shCoefficients[ 1 ] * 2.0 * 0.511664 * y;
	result += shCoefficients[ 2 ] * 2.0 * 0.511664 * z;
	result += shCoefficients[ 3 ] * 2.0 * 0.511664 * x;
	result += shCoefficients[ 4 ] * 2.0 * 0.429043 * x * y;
	result += shCoefficients[ 5 ] * 2.0 * 0.429043 * y * z;
	result += shCoefficients[ 6 ] * ( 0.743125 * z * z - 0.247708 );
	result += shCoefficients[ 7 ] * 2.0 * 0.429043 * x * z;
	result += shCoefficients[ 8 ] * 0.429043 * ( x * x - y * y );
	return result;
}
vec3 getLightProbeIrradiance( const in vec3 lightProbe[ 9 ], const in vec3 normal ) {
	vec3 worldNormal = inverseTransformDirection( normal, viewMatrix );
	vec3 irradiance = shGetIrradianceAt( worldNormal, lightProbe );
	return irradiance;
}
vec3 getAmbientLightIrradiance( const in vec3 ambientLightColor ) {
	vec3 irradiance = ambientLightColor;
	return irradiance;
}
float getDistanceAttenuation( const in float lightDistance, const in float cutoffDistance, const in float decayExponent ) {
	#if defined ( LEGACY_LIGHTS )
		if ( cutoffDistance > 0.0 && decayExponent > 0.0 ) {
			return pow( saturate( - lightDistance / cutoffDistance + 1.0 ), decayExponent );
		}
		return 1.0;
	#else
		float distanceFalloff = 1.0 / max( pow( lightDistance, decayExponent ), 0.01 );
		if ( cutoffDistance > 0.0 ) {
			distanceFalloff *= pow2( saturate( 1.0 - pow4( lightDistance / cutoffDistance ) ) );
		}
		return distanceFalloff;
	#endif
}
float getSpotAttenuation( const in float coneCosine, const in float penumbraCosine, const in float angleCosine ) {
	return smoothstep( coneCosine, penumbraCosine, angleCosine );
}
#if NUM_DIR_LIGHTS > 0
	struct DirectionalLight {
		vec3 direction;
		vec3 color;
	};
	uniform DirectionalLight directionalLights[ NUM_DIR_LIGHTS ];
	void getDirectionalLightInfo( const in DirectionalLight directionalLight, const in GeometricContext geometry, out IncidentLight light ) {
		light.color = directionalLight.color;
		light.direction = directionalLight.direction;
		light.visible = true;
	}
#endif
#if NUM_POINT_LIGHTS > 0
	struct PointLight {
		vec3 position;
		vec3 color;
		float distance;
		float decay;
	};
	uniform PointLight pointLights[ NUM_POINT_LIGHTS ];
	void getPointLightInfo( const in PointLight pointLight, const in GeometricContext geometry, out IncidentLight light ) {
		vec3 lVector = pointLight.position - geometry.position;
		light.direction = normalize( lVector );
		float lightDistance = length( lVector );
		light.color = pointLight.color;
		light.color *= getDistanceAttenuation( lightDistance, pointLight.distance, pointLight.decay );
		light.visible = ( light.color != vec3( 0.0 ) );
	}
#endif
#if NUM_SPOT_LIGHTS > 0
	struct SpotLight {
		vec3 position;
		vec3 direction;
		vec3 color;
		float distance;
		float decay;
		float coneCos;
		float penumbraCos;
	};
	uniform SpotLight spotLights[ NUM_SPOT_LIGHTS ];
	void getSpotLightInfo( const in SpotLight spotLight, const in GeometricContext geometry, out IncidentLight light ) {
		vec3 lVector = spotLight.position - geometry.position;
		light.direction = normalize( lVector );
		float angleCos = dot( light.direction, spotLight.direction );
		float spotAttenuation = getSpotAttenuation( spotLight.coneCos, spotLight.penumbraCos, angleCos );
		if ( spotAttenuation > 0.0 ) {
			float lightDistance = length( lVector );
			light.color = spotLight.color * spotAttenuation;
			light.color *= getDistanceAttenuation( lightDistance, spotLight.distance, spotLight.decay );
			light.visible = ( light.color != vec3( 0.0 ) );
		} else {
			light.color = vec3( 0.0 );
			light.visible = false;
		}
	}
#endif
#if NUM_RECT_AREA_LIGHTS > 0
	struct RectAreaLight {
		vec3 color;
		vec3 position;
		vec3 halfWidth;
		vec3 halfHeight;
	};
	uniform sampler2D ltc_1;	uniform sampler2D ltc_2;
	uniform RectAreaLight rectAreaLights[ NUM_RECT_AREA_LIGHTS ];
#endif
#if NUM_HEMI_LIGHTS > 0
	struct HemisphereLight {
		vec3 direction;
		vec3 skyColor;
		vec3 groundColor;
	};
	uniform HemisphereLight hemisphereLights[ NUM_HEMI_LIGHTS ];
	vec3 getHemisphereLightIrradiance( const in HemisphereLight hemiLight, const in vec3 normal ) {
		float dotNL = dot( normal, hemiLight.direction );
		float hemiDiffuseWeight = 0.5 * dotNL + 0.5;
		vec3 irradiance = mix( hemiLight.groundColor, hemiLight.skyColor, hemiDiffuseWeight );
		return irradiance;
	}
#endif`, wh = `#ifdef USE_ENVMAP
	vec3 getIBLIrradiance( const in vec3 normal ) {
		#ifdef ENVMAP_TYPE_CUBE_UV
			vec3 worldNormal = inverseTransformDirection( normal, viewMatrix );
			vec4 envMapColor = textureCubeUV( envMap, worldNormal, 1.0 );
			return PI * envMapColor.rgb * envMapIntensity;
		#else
			return vec3( 0.0 );
		#endif
	}
	vec3 getIBLRadiance( const in vec3 viewDir, const in vec3 normal, const in float roughness ) {
		#ifdef ENVMAP_TYPE_CUBE_UV
			vec3 reflectVec = reflect( - viewDir, normal );
			reflectVec = normalize( mix( reflectVec, normal, roughness * roughness) );
			reflectVec = inverseTransformDirection( reflectVec, viewMatrix );
			vec4 envMapColor = textureCubeUV( envMap, reflectVec, roughness );
			return envMapColor.rgb * envMapIntensity;
		#else
			return vec3( 0.0 );
		#endif
	}
	#ifdef USE_ANISOTROPY
		vec3 getIBLAnisotropyRadiance( const in vec3 viewDir, const in vec3 normal, const in float roughness, const in vec3 bitangent, const in float anisotropy ) {
			#ifdef ENVMAP_TYPE_CUBE_UV
				vec3 bentNormal = cross( bitangent, viewDir );
				bentNormal = normalize( cross( bentNormal, bitangent ) );
				bentNormal = normalize( mix( bentNormal, normal, pow2( pow2( 1.0 - anisotropy * ( 1.0 - roughness ) ) ) ) );
				return getIBLRadiance( viewDir, bentNormal, roughness );
			#else
				return vec3( 0.0 );
			#endif
		}
	#endif
#endif`, Rh = `ToonMaterial material;
material.diffuseColor = diffuseColor.rgb;`, Ch = `varying vec3 vViewPosition;
struct ToonMaterial {
	vec3 diffuseColor;
};
void RE_Direct_Toon( const in IncidentLight directLight, const in GeometricContext geometry, const in ToonMaterial material, inout ReflectedLight reflectedLight ) {
	vec3 irradiance = getGradientIrradiance( geometry.normal, directLight.direction ) * directLight.color;
	reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
void RE_IndirectDiffuse_Toon( const in vec3 irradiance, const in GeometricContext geometry, const in ToonMaterial material, inout ReflectedLight reflectedLight ) {
	reflectedLight.indirectDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
#define RE_Direct				RE_Direct_Toon
#define RE_IndirectDiffuse		RE_IndirectDiffuse_Toon`, Ph = `BlinnPhongMaterial material;
material.diffuseColor = diffuseColor.rgb;
material.specularColor = specular;
material.specularShininess = shininess;
material.specularStrength = specularStrength;`, Lh = `varying vec3 vViewPosition;
struct BlinnPhongMaterial {
	vec3 diffuseColor;
	vec3 specularColor;
	float specularShininess;
	float specularStrength;
};
void RE_Direct_BlinnPhong( const in IncidentLight directLight, const in GeometricContext geometry, const in BlinnPhongMaterial material, inout ReflectedLight reflectedLight ) {
	float dotNL = saturate( dot( geometry.normal, directLight.direction ) );
	vec3 irradiance = dotNL * directLight.color;
	reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
	reflectedLight.directSpecular += irradiance * BRDF_BlinnPhong( directLight.direction, geometry.viewDir, geometry.normal, material.specularColor, material.specularShininess ) * material.specularStrength;
}
void RE_IndirectDiffuse_BlinnPhong( const in vec3 irradiance, const in GeometricContext geometry, const in BlinnPhongMaterial material, inout ReflectedLight reflectedLight ) {
	reflectedLight.indirectDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
#define RE_Direct				RE_Direct_BlinnPhong
#define RE_IndirectDiffuse		RE_IndirectDiffuse_BlinnPhong`, Uh = `PhysicalMaterial material;
material.diffuseColor = diffuseColor.rgb * ( 1.0 - metalnessFactor );
vec3 dxy = max( abs( dFdx( geometryNormal ) ), abs( dFdy( geometryNormal ) ) );
float geometryRoughness = max( max( dxy.x, dxy.y ), dxy.z );
material.roughness = max( roughnessFactor, 0.0525 );material.roughness += geometryRoughness;
material.roughness = min( material.roughness, 1.0 );
#ifdef IOR
	material.ior = ior;
	#ifdef USE_SPECULAR
		float specularIntensityFactor = specularIntensity;
		vec3 specularColorFactor = specularColor;
		#ifdef USE_SPECULAR_COLORMAP
			specularColorFactor *= texture2D( specularColorMap, vSpecularColorMapUv ).rgb;
		#endif
		#ifdef USE_SPECULAR_INTENSITYMAP
			specularIntensityFactor *= texture2D( specularIntensityMap, vSpecularIntensityMapUv ).a;
		#endif
		material.specularF90 = mix( specularIntensityFactor, 1.0, metalnessFactor );
	#else
		float specularIntensityFactor = 1.0;
		vec3 specularColorFactor = vec3( 1.0 );
		material.specularF90 = 1.0;
	#endif
	material.specularColor = mix( min( pow2( ( material.ior - 1.0 ) / ( material.ior + 1.0 ) ) * specularColorFactor, vec3( 1.0 ) ) * specularIntensityFactor, diffuseColor.rgb, metalnessFactor );
#else
	material.specularColor = mix( vec3( 0.04 ), diffuseColor.rgb, metalnessFactor );
	material.specularF90 = 1.0;
#endif
#ifdef USE_CLEARCOAT
	material.clearcoat = clearcoat;
	material.clearcoatRoughness = clearcoatRoughness;
	material.clearcoatF0 = vec3( 0.04 );
	material.clearcoatF90 = 1.0;
	#ifdef USE_CLEARCOATMAP
		material.clearcoat *= texture2D( clearcoatMap, vClearcoatMapUv ).x;
	#endif
	#ifdef USE_CLEARCOAT_ROUGHNESSMAP
		material.clearcoatRoughness *= texture2D( clearcoatRoughnessMap, vClearcoatRoughnessMapUv ).y;
	#endif
	material.clearcoat = saturate( material.clearcoat );	material.clearcoatRoughness = max( material.clearcoatRoughness, 0.0525 );
	material.clearcoatRoughness += geometryRoughness;
	material.clearcoatRoughness = min( material.clearcoatRoughness, 1.0 );
#endif
#ifdef USE_IRIDESCENCE
	material.iridescence = iridescence;
	material.iridescenceIOR = iridescenceIOR;
	#ifdef USE_IRIDESCENCEMAP
		material.iridescence *= texture2D( iridescenceMap, vIridescenceMapUv ).r;
	#endif
	#ifdef USE_IRIDESCENCE_THICKNESSMAP
		material.iridescenceThickness = (iridescenceThicknessMaximum - iridescenceThicknessMinimum) * texture2D( iridescenceThicknessMap, vIridescenceThicknessMapUv ).g + iridescenceThicknessMinimum;
	#else
		material.iridescenceThickness = iridescenceThicknessMaximum;
	#endif
#endif
#ifdef USE_SHEEN
	material.sheenColor = sheenColor;
	#ifdef USE_SHEEN_COLORMAP
		material.sheenColor *= texture2D( sheenColorMap, vSheenColorMapUv ).rgb;
	#endif
	material.sheenRoughness = clamp( sheenRoughness, 0.07, 1.0 );
	#ifdef USE_SHEEN_ROUGHNESSMAP
		material.sheenRoughness *= texture2D( sheenRoughnessMap, vSheenRoughnessMapUv ).a;
	#endif
#endif
#ifdef USE_ANISOTROPY
	#ifdef USE_ANISOTROPYMAP
		mat2 anisotropyMat = mat2( anisotropyVector.x, anisotropyVector.y, - anisotropyVector.y, anisotropyVector.x );
		vec3 anisotropyPolar = texture2D( anisotropyMap, vAnisotropyMapUv ).rgb;
		vec2 anisotropyV = anisotropyMat * normalize( 2.0 * anisotropyPolar.rg - vec2( 1.0 ) ) * anisotropyPolar.b;
	#else
		vec2 anisotropyV = anisotropyVector;
	#endif
	material.anisotropy = length( anisotropyV );
	anisotropyV /= material.anisotropy;
	material.anisotropy = saturate( material.anisotropy );
	material.alphaT = mix( pow2( material.roughness ), 1.0, pow2( material.anisotropy ) );
	material.anisotropyT = tbn[ 0 ] * anisotropyV.x - tbn[ 1 ] * anisotropyV.y;
	material.anisotropyB = tbn[ 1 ] * anisotropyV.x + tbn[ 0 ] * anisotropyV.y;
#endif`, Dh = `struct PhysicalMaterial {
	vec3 diffuseColor;
	float roughness;
	vec3 specularColor;
	float specularF90;
	#ifdef USE_CLEARCOAT
		float clearcoat;
		float clearcoatRoughness;
		vec3 clearcoatF0;
		float clearcoatF90;
	#endif
	#ifdef USE_IRIDESCENCE
		float iridescence;
		float iridescenceIOR;
		float iridescenceThickness;
		vec3 iridescenceFresnel;
		vec3 iridescenceF0;
	#endif
	#ifdef USE_SHEEN
		vec3 sheenColor;
		float sheenRoughness;
	#endif
	#ifdef IOR
		float ior;
	#endif
	#ifdef USE_TRANSMISSION
		float transmission;
		float transmissionAlpha;
		float thickness;
		float attenuationDistance;
		vec3 attenuationColor;
	#endif
	#ifdef USE_ANISOTROPY
		float anisotropy;
		float alphaT;
		vec3 anisotropyT;
		vec3 anisotropyB;
	#endif
};
vec3 clearcoatSpecular = vec3( 0.0 );
vec3 sheenSpecular = vec3( 0.0 );
vec3 Schlick_to_F0( const in vec3 f, const in float f90, const in float dotVH ) {
    float x = clamp( 1.0 - dotVH, 0.0, 1.0 );
    float x2 = x * x;
    float x5 = clamp( x * x2 * x2, 0.0, 0.9999 );
    return ( f - vec3( f90 ) * x5 ) / ( 1.0 - x5 );
}
float V_GGX_SmithCorrelated( const in float alpha, const in float dotNL, const in float dotNV ) {
	float a2 = pow2( alpha );
	float gv = dotNL * sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNV ) );
	float gl = dotNV * sqrt( a2 + ( 1.0 - a2 ) * pow2( dotNL ) );
	return 0.5 / max( gv + gl, EPSILON );
}
float D_GGX( const in float alpha, const in float dotNH ) {
	float a2 = pow2( alpha );
	float denom = pow2( dotNH ) * ( a2 - 1.0 ) + 1.0;
	return RECIPROCAL_PI * a2 / pow2( denom );
}
#ifdef USE_ANISOTROPY
	float V_GGX_SmithCorrelated_Anisotropic( const in float alphaT, const in float alphaB, const in float dotTV, const in float dotBV, const in float dotTL, const in float dotBL, const in float dotNV, const in float dotNL ) {
		float gv = dotNL * length( vec3( alphaT * dotTV, alphaB * dotBV, dotNV ) );
		float gl = dotNV * length( vec3( alphaT * dotTL, alphaB * dotBL, dotNL ) );
		float v = 0.5 / ( gv + gl );
		return saturate(v);
	}
	float D_GGX_Anisotropic( const in float alphaT, const in float alphaB, const in float dotNH, const in float dotTH, const in float dotBH ) {
		float a2 = alphaT * alphaB;
		highp vec3 v = vec3( alphaB * dotTH, alphaT * dotBH, a2 * dotNH );
		highp float v2 = dot( v, v );
		float w2 = a2 / v2;
		return RECIPROCAL_PI * a2 * pow2 ( w2 );
	}
#endif
#ifdef USE_CLEARCOAT
	vec3 BRDF_GGX_Clearcoat( const in vec3 lightDir, const in vec3 viewDir, const in vec3 normal, const in PhysicalMaterial material) {
		vec3 f0 = material.clearcoatF0;
		float f90 = material.clearcoatF90;
		float roughness = material.clearcoatRoughness;
		float alpha = pow2( roughness );
		vec3 halfDir = normalize( lightDir + viewDir );
		float dotNL = saturate( dot( normal, lightDir ) );
		float dotNV = saturate( dot( normal, viewDir ) );
		float dotNH = saturate( dot( normal, halfDir ) );
		float dotVH = saturate( dot( viewDir, halfDir ) );
		vec3 F = F_Schlick( f0, f90, dotVH );
		float V = V_GGX_SmithCorrelated( alpha, dotNL, dotNV );
		float D = D_GGX( alpha, dotNH );
		return F * ( V * D );
	}
#endif
vec3 BRDF_GGX( const in vec3 lightDir, const in vec3 viewDir, const in vec3 normal, const in PhysicalMaterial material ) {
	vec3 f0 = material.specularColor;
	float f90 = material.specularF90;
	float roughness = material.roughness;
	float alpha = pow2( roughness );
	vec3 halfDir = normalize( lightDir + viewDir );
	float dotNL = saturate( dot( normal, lightDir ) );
	float dotNV = saturate( dot( normal, viewDir ) );
	float dotNH = saturate( dot( normal, halfDir ) );
	float dotVH = saturate( dot( viewDir, halfDir ) );
	vec3 F = F_Schlick( f0, f90, dotVH );
	#ifdef USE_IRIDESCENCE
		F = mix( F, material.iridescenceFresnel, material.iridescence );
	#endif
	#ifdef USE_ANISOTROPY
		float dotTL = dot( material.anisotropyT, lightDir );
		float dotTV = dot( material.anisotropyT, viewDir );
		float dotTH = dot( material.anisotropyT, halfDir );
		float dotBL = dot( material.anisotropyB, lightDir );
		float dotBV = dot( material.anisotropyB, viewDir );
		float dotBH = dot( material.anisotropyB, halfDir );
		float V = V_GGX_SmithCorrelated_Anisotropic( material.alphaT, alpha, dotTV, dotBV, dotTL, dotBL, dotNV, dotNL );
		float D = D_GGX_Anisotropic( material.alphaT, alpha, dotNH, dotTH, dotBH );
	#else
		float V = V_GGX_SmithCorrelated( alpha, dotNL, dotNV );
		float D = D_GGX( alpha, dotNH );
	#endif
	return F * ( V * D );
}
vec2 LTC_Uv( const in vec3 N, const in vec3 V, const in float roughness ) {
	const float LUT_SIZE = 64.0;
	const float LUT_SCALE = ( LUT_SIZE - 1.0 ) / LUT_SIZE;
	const float LUT_BIAS = 0.5 / LUT_SIZE;
	float dotNV = saturate( dot( N, V ) );
	vec2 uv = vec2( roughness, sqrt( 1.0 - dotNV ) );
	uv = uv * LUT_SCALE + LUT_BIAS;
	return uv;
}
float LTC_ClippedSphereFormFactor( const in vec3 f ) {
	float l = length( f );
	return max( ( l * l + f.z ) / ( l + 1.0 ), 0.0 );
}
vec3 LTC_EdgeVectorFormFactor( const in vec3 v1, const in vec3 v2 ) {
	float x = dot( v1, v2 );
	float y = abs( x );
	float a = 0.8543985 + ( 0.4965155 + 0.0145206 * y ) * y;
	float b = 3.4175940 + ( 4.1616724 + y ) * y;
	float v = a / b;
	float theta_sintheta = ( x > 0.0 ) ? v : 0.5 * inversesqrt( max( 1.0 - x * x, 1e-7 ) ) - v;
	return cross( v1, v2 ) * theta_sintheta;
}
vec3 LTC_Evaluate( const in vec3 N, const in vec3 V, const in vec3 P, const in mat3 mInv, const in vec3 rectCoords[ 4 ] ) {
	vec3 v1 = rectCoords[ 1 ] - rectCoords[ 0 ];
	vec3 v2 = rectCoords[ 3 ] - rectCoords[ 0 ];
	vec3 lightNormal = cross( v1, v2 );
	if( dot( lightNormal, P - rectCoords[ 0 ] ) < 0.0 ) return vec3( 0.0 );
	vec3 T1, T2;
	T1 = normalize( V - N * dot( V, N ) );
	T2 = - cross( N, T1 );
	mat3 mat = mInv * transposeMat3( mat3( T1, T2, N ) );
	vec3 coords[ 4 ];
	coords[ 0 ] = mat * ( rectCoords[ 0 ] - P );
	coords[ 1 ] = mat * ( rectCoords[ 1 ] - P );
	coords[ 2 ] = mat * ( rectCoords[ 2 ] - P );
	coords[ 3 ] = mat * ( rectCoords[ 3 ] - P );
	coords[ 0 ] = normalize( coords[ 0 ] );
	coords[ 1 ] = normalize( coords[ 1 ] );
	coords[ 2 ] = normalize( coords[ 2 ] );
	coords[ 3 ] = normalize( coords[ 3 ] );
	vec3 vectorFormFactor = vec3( 0.0 );
	vectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 0 ], coords[ 1 ] );
	vectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 1 ], coords[ 2 ] );
	vectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 2 ], coords[ 3 ] );
	vectorFormFactor += LTC_EdgeVectorFormFactor( coords[ 3 ], coords[ 0 ] );
	float result = LTC_ClippedSphereFormFactor( vectorFormFactor );
	return vec3( result );
}
#if defined( USE_SHEEN )
float D_Charlie( float roughness, float dotNH ) {
	float alpha = pow2( roughness );
	float invAlpha = 1.0 / alpha;
	float cos2h = dotNH * dotNH;
	float sin2h = max( 1.0 - cos2h, 0.0078125 );
	return ( 2.0 + invAlpha ) * pow( sin2h, invAlpha * 0.5 ) / ( 2.0 * PI );
}
float V_Neubelt( float dotNV, float dotNL ) {
	return saturate( 1.0 / ( 4.0 * ( dotNL + dotNV - dotNL * dotNV ) ) );
}
vec3 BRDF_Sheen( const in vec3 lightDir, const in vec3 viewDir, const in vec3 normal, vec3 sheenColor, const in float sheenRoughness ) {
	vec3 halfDir = normalize( lightDir + viewDir );
	float dotNL = saturate( dot( normal, lightDir ) );
	float dotNV = saturate( dot( normal, viewDir ) );
	float dotNH = saturate( dot( normal, halfDir ) );
	float D = D_Charlie( sheenRoughness, dotNH );
	float V = V_Neubelt( dotNV, dotNL );
	return sheenColor * ( D * V );
}
#endif
float IBLSheenBRDF( const in vec3 normal, const in vec3 viewDir, const in float roughness ) {
	float dotNV = saturate( dot( normal, viewDir ) );
	float r2 = roughness * roughness;
	float a = roughness < 0.25 ? -339.2 * r2 + 161.4 * roughness - 25.9 : -8.48 * r2 + 14.3 * roughness - 9.95;
	float b = roughness < 0.25 ? 44.0 * r2 - 23.7 * roughness + 3.26 : 1.97 * r2 - 3.27 * roughness + 0.72;
	float DG = exp( a * dotNV + b ) + ( roughness < 0.25 ? 0.0 : 0.1 * ( roughness - 0.25 ) );
	return saturate( DG * RECIPROCAL_PI );
}
vec2 DFGApprox( const in vec3 normal, const in vec3 viewDir, const in float roughness ) {
	float dotNV = saturate( dot( normal, viewDir ) );
	const vec4 c0 = vec4( - 1, - 0.0275, - 0.572, 0.022 );
	const vec4 c1 = vec4( 1, 0.0425, 1.04, - 0.04 );
	vec4 r = roughness * c0 + c1;
	float a004 = min( r.x * r.x, exp2( - 9.28 * dotNV ) ) * r.x + r.y;
	vec2 fab = vec2( - 1.04, 1.04 ) * a004 + r.zw;
	return fab;
}
vec3 EnvironmentBRDF( const in vec3 normal, const in vec3 viewDir, const in vec3 specularColor, const in float specularF90, const in float roughness ) {
	vec2 fab = DFGApprox( normal, viewDir, roughness );
	return specularColor * fab.x + specularF90 * fab.y;
}
#ifdef USE_IRIDESCENCE
void computeMultiscatteringIridescence( const in vec3 normal, const in vec3 viewDir, const in vec3 specularColor, const in float specularF90, const in float iridescence, const in vec3 iridescenceF0, const in float roughness, inout vec3 singleScatter, inout vec3 multiScatter ) {
#else
void computeMultiscattering( const in vec3 normal, const in vec3 viewDir, const in vec3 specularColor, const in float specularF90, const in float roughness, inout vec3 singleScatter, inout vec3 multiScatter ) {
#endif
	vec2 fab = DFGApprox( normal, viewDir, roughness );
	#ifdef USE_IRIDESCENCE
		vec3 Fr = mix( specularColor, iridescenceF0, iridescence );
	#else
		vec3 Fr = specularColor;
	#endif
	vec3 FssEss = Fr * fab.x + specularF90 * fab.y;
	float Ess = fab.x + fab.y;
	float Ems = 1.0 - Ess;
	vec3 Favg = Fr + ( 1.0 - Fr ) * 0.047619;	vec3 Fms = FssEss * Favg / ( 1.0 - Ems * Favg );
	singleScatter += FssEss;
	multiScatter += Fms * Ems;
}
#if NUM_RECT_AREA_LIGHTS > 0
	void RE_Direct_RectArea_Physical( const in RectAreaLight rectAreaLight, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {
		vec3 normal = geometry.normal;
		vec3 viewDir = geometry.viewDir;
		vec3 position = geometry.position;
		vec3 lightPos = rectAreaLight.position;
		vec3 halfWidth = rectAreaLight.halfWidth;
		vec3 halfHeight = rectAreaLight.halfHeight;
		vec3 lightColor = rectAreaLight.color;
		float roughness = material.roughness;
		vec3 rectCoords[ 4 ];
		rectCoords[ 0 ] = lightPos + halfWidth - halfHeight;		rectCoords[ 1 ] = lightPos - halfWidth - halfHeight;
		rectCoords[ 2 ] = lightPos - halfWidth + halfHeight;
		rectCoords[ 3 ] = lightPos + halfWidth + halfHeight;
		vec2 uv = LTC_Uv( normal, viewDir, roughness );
		vec4 t1 = texture2D( ltc_1, uv );
		vec4 t2 = texture2D( ltc_2, uv );
		mat3 mInv = mat3(
			vec3( t1.x, 0, t1.y ),
			vec3(    0, 1,    0 ),
			vec3( t1.z, 0, t1.w )
		);
		vec3 fresnel = ( material.specularColor * t2.x + ( vec3( 1.0 ) - material.specularColor ) * t2.y );
		reflectedLight.directSpecular += lightColor * fresnel * LTC_Evaluate( normal, viewDir, position, mInv, rectCoords );
		reflectedLight.directDiffuse += lightColor * material.diffuseColor * LTC_Evaluate( normal, viewDir, position, mat3( 1.0 ), rectCoords );
	}
#endif
void RE_Direct_Physical( const in IncidentLight directLight, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {
	float dotNL = saturate( dot( geometry.normal, directLight.direction ) );
	vec3 irradiance = dotNL * directLight.color;
	#ifdef USE_CLEARCOAT
		float dotNLcc = saturate( dot( geometry.clearcoatNormal, directLight.direction ) );
		vec3 ccIrradiance = dotNLcc * directLight.color;
		clearcoatSpecular += ccIrradiance * BRDF_GGX_Clearcoat( directLight.direction, geometry.viewDir, geometry.clearcoatNormal, material );
	#endif
	#ifdef USE_SHEEN
		sheenSpecular += irradiance * BRDF_Sheen( directLight.direction, geometry.viewDir, geometry.normal, material.sheenColor, material.sheenRoughness );
	#endif
	reflectedLight.directSpecular += irradiance * BRDF_GGX( directLight.direction, geometry.viewDir, geometry.normal, material );
	reflectedLight.directDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
void RE_IndirectDiffuse_Physical( const in vec3 irradiance, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight ) {
	reflectedLight.indirectDiffuse += irradiance * BRDF_Lambert( material.diffuseColor );
}
void RE_IndirectSpecular_Physical( const in vec3 radiance, const in vec3 irradiance, const in vec3 clearcoatRadiance, const in GeometricContext geometry, const in PhysicalMaterial material, inout ReflectedLight reflectedLight) {
	#ifdef USE_CLEARCOAT
		clearcoatSpecular += clearcoatRadiance * EnvironmentBRDF( geometry.clearcoatNormal, geometry.viewDir, material.clearcoatF0, material.clearcoatF90, material.clearcoatRoughness );
	#endif
	#ifdef USE_SHEEN
		sheenSpecular += irradiance * material.sheenColor * IBLSheenBRDF( geometry.normal, geometry.viewDir, material.sheenRoughness );
	#endif
	vec3 singleScattering = vec3( 0.0 );
	vec3 multiScattering = vec3( 0.0 );
	vec3 cosineWeightedIrradiance = irradiance * RECIPROCAL_PI;
	#ifdef USE_IRIDESCENCE
		computeMultiscatteringIridescence( geometry.normal, geometry.viewDir, material.specularColor, material.specularF90, material.iridescence, material.iridescenceFresnel, material.roughness, singleScattering, multiScattering );
	#else
		computeMultiscattering( geometry.normal, geometry.viewDir, material.specularColor, material.specularF90, material.roughness, singleScattering, multiScattering );
	#endif
	vec3 totalScattering = singleScattering + multiScattering;
	vec3 diffuse = material.diffuseColor * ( 1.0 - max( max( totalScattering.r, totalScattering.g ), totalScattering.b ) );
	reflectedLight.indirectSpecular += radiance * singleScattering;
	reflectedLight.indirectSpecular += multiScattering * cosineWeightedIrradiance;
	reflectedLight.indirectDiffuse += diffuse * cosineWeightedIrradiance;
}
#define RE_Direct				RE_Direct_Physical
#define RE_Direct_RectArea		RE_Direct_RectArea_Physical
#define RE_IndirectDiffuse		RE_IndirectDiffuse_Physical
#define RE_IndirectSpecular		RE_IndirectSpecular_Physical
float computeSpecularOcclusion( const in float dotNV, const in float ambientOcclusion, const in float roughness ) {
	return saturate( pow( dotNV + ambientOcclusion, exp2( - 16.0 * roughness - 1.0 ) ) - 1.0 + ambientOcclusion );
}`, Ih = `
GeometricContext geometry;
geometry.position = - vViewPosition;
geometry.normal = normal;
geometry.viewDir = ( isOrthographic ) ? vec3( 0, 0, 1 ) : normalize( vViewPosition );
#ifdef USE_CLEARCOAT
	geometry.clearcoatNormal = clearcoatNormal;
#endif
#ifdef USE_IRIDESCENCE
	float dotNVi = saturate( dot( normal, geometry.viewDir ) );
	if ( material.iridescenceThickness == 0.0 ) {
		material.iridescence = 0.0;
	} else {
		material.iridescence = saturate( material.iridescence );
	}
	if ( material.iridescence > 0.0 ) {
		material.iridescenceFresnel = evalIridescence( 1.0, material.iridescenceIOR, dotNVi, material.iridescenceThickness, material.specularColor );
		material.iridescenceF0 = Schlick_to_F0( material.iridescenceFresnel, 1.0, dotNVi );
	}
#endif
IncidentLight directLight;
#if ( NUM_POINT_LIGHTS > 0 ) && defined( RE_Direct )
	PointLight pointLight;
	#if defined( USE_SHADOWMAP ) && NUM_POINT_LIGHT_SHADOWS > 0
	PointLightShadow pointLightShadow;
	#endif
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_POINT_LIGHTS; i ++ ) {
		pointLight = pointLights[ i ];
		getPointLightInfo( pointLight, geometry, directLight );
		#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_POINT_LIGHT_SHADOWS )
		pointLightShadow = pointLightShadows[ i ];
		directLight.color *= ( directLight.visible && receiveShadow ) ? getPointShadow( pointShadowMap[ i ], pointLightShadow.shadowMapSize, pointLightShadow.shadowBias, pointLightShadow.shadowRadius, vPointShadowCoord[ i ], pointLightShadow.shadowCameraNear, pointLightShadow.shadowCameraFar ) : 1.0;
		#endif
		RE_Direct( directLight, geometry, material, reflectedLight );
	}
	#pragma unroll_loop_end
#endif
#if ( NUM_SPOT_LIGHTS > 0 ) && defined( RE_Direct )
	SpotLight spotLight;
	vec4 spotColor;
	vec3 spotLightCoord;
	bool inSpotLightMap;
	#if defined( USE_SHADOWMAP ) && NUM_SPOT_LIGHT_SHADOWS > 0
	SpotLightShadow spotLightShadow;
	#endif
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_SPOT_LIGHTS; i ++ ) {
		spotLight = spotLights[ i ];
		getSpotLightInfo( spotLight, geometry, directLight );
		#if ( UNROLLED_LOOP_INDEX < NUM_SPOT_LIGHT_SHADOWS_WITH_MAPS )
		#define SPOT_LIGHT_MAP_INDEX UNROLLED_LOOP_INDEX
		#elif ( UNROLLED_LOOP_INDEX < NUM_SPOT_LIGHT_SHADOWS )
		#define SPOT_LIGHT_MAP_INDEX NUM_SPOT_LIGHT_MAPS
		#else
		#define SPOT_LIGHT_MAP_INDEX ( UNROLLED_LOOP_INDEX - NUM_SPOT_LIGHT_SHADOWS + NUM_SPOT_LIGHT_SHADOWS_WITH_MAPS )
		#endif
		#if ( SPOT_LIGHT_MAP_INDEX < NUM_SPOT_LIGHT_MAPS )
			spotLightCoord = vSpotLightCoord[ i ].xyz / vSpotLightCoord[ i ].w;
			inSpotLightMap = all( lessThan( abs( spotLightCoord * 2. - 1. ), vec3( 1.0 ) ) );
			spotColor = texture2D( spotLightMap[ SPOT_LIGHT_MAP_INDEX ], spotLightCoord.xy );
			directLight.color = inSpotLightMap ? directLight.color * spotColor.rgb : directLight.color;
		#endif
		#undef SPOT_LIGHT_MAP_INDEX
		#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_SPOT_LIGHT_SHADOWS )
		spotLightShadow = spotLightShadows[ i ];
		directLight.color *= ( directLight.visible && receiveShadow ) ? getShadow( spotShadowMap[ i ], spotLightShadow.shadowMapSize, spotLightShadow.shadowBias, spotLightShadow.shadowRadius, vSpotLightCoord[ i ] ) : 1.0;
		#endif
		RE_Direct( directLight, geometry, material, reflectedLight );
	}
	#pragma unroll_loop_end
#endif
#if ( NUM_DIR_LIGHTS > 0 ) && defined( RE_Direct )
	DirectionalLight directionalLight;
	#if defined( USE_SHADOWMAP ) && NUM_DIR_LIGHT_SHADOWS > 0
	DirectionalLightShadow directionalLightShadow;
	#endif
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_DIR_LIGHTS; i ++ ) {
		directionalLight = directionalLights[ i ];
		getDirectionalLightInfo( directionalLight, geometry, directLight );
		#if defined( USE_SHADOWMAP ) && ( UNROLLED_LOOP_INDEX < NUM_DIR_LIGHT_SHADOWS )
		directionalLightShadow = directionalLightShadows[ i ];
		directLight.color *= ( directLight.visible && receiveShadow ) ? getShadow( directionalShadowMap[ i ], directionalLightShadow.shadowMapSize, directionalLightShadow.shadowBias, directionalLightShadow.shadowRadius, vDirectionalShadowCoord[ i ] ) : 1.0;
		#endif
		RE_Direct( directLight, geometry, material, reflectedLight );
	}
	#pragma unroll_loop_end
#endif
#if ( NUM_RECT_AREA_LIGHTS > 0 ) && defined( RE_Direct_RectArea )
	RectAreaLight rectAreaLight;
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_RECT_AREA_LIGHTS; i ++ ) {
		rectAreaLight = rectAreaLights[ i ];
		RE_Direct_RectArea( rectAreaLight, geometry, material, reflectedLight );
	}
	#pragma unroll_loop_end
#endif
#if defined( RE_IndirectDiffuse )
	vec3 iblIrradiance = vec3( 0.0 );
	vec3 irradiance = getAmbientLightIrradiance( ambientLightColor );
	irradiance += getLightProbeIrradiance( lightProbe, geometry.normal );
	#if ( NUM_HEMI_LIGHTS > 0 )
		#pragma unroll_loop_start
		for ( int i = 0; i < NUM_HEMI_LIGHTS; i ++ ) {
			irradiance += getHemisphereLightIrradiance( hemisphereLights[ i ], geometry.normal );
		}
		#pragma unroll_loop_end
	#endif
#endif
#if defined( RE_IndirectSpecular )
	vec3 radiance = vec3( 0.0 );
	vec3 clearcoatRadiance = vec3( 0.0 );
#endif`, Nh = `#if defined( RE_IndirectDiffuse )
	#ifdef USE_LIGHTMAP
		vec4 lightMapTexel = texture2D( lightMap, vLightMapUv );
		vec3 lightMapIrradiance = lightMapTexel.rgb * lightMapIntensity;
		irradiance += lightMapIrradiance;
	#endif
	#if defined( USE_ENVMAP ) && defined( STANDARD ) && defined( ENVMAP_TYPE_CUBE_UV )
		iblIrradiance += getIBLIrradiance( geometry.normal );
	#endif
#endif
#if defined( USE_ENVMAP ) && defined( RE_IndirectSpecular )
	#ifdef USE_ANISOTROPY
		radiance += getIBLAnisotropyRadiance( geometry.viewDir, geometry.normal, material.roughness, material.anisotropyB, material.anisotropy );
	#else
		radiance += getIBLRadiance( geometry.viewDir, geometry.normal, material.roughness );
	#endif
	#ifdef USE_CLEARCOAT
		clearcoatRadiance += getIBLRadiance( geometry.viewDir, geometry.clearcoatNormal, material.clearcoatRoughness );
	#endif
#endif`, Oh = `#if defined( RE_IndirectDiffuse )
	RE_IndirectDiffuse( irradiance, geometry, material, reflectedLight );
#endif
#if defined( RE_IndirectSpecular )
	RE_IndirectSpecular( radiance, iblIrradiance, clearcoatRadiance, geometry, material, reflectedLight );
#endif`, Fh = `#if defined( USE_LOGDEPTHBUF ) && defined( USE_LOGDEPTHBUF_EXT )
	gl_FragDepthEXT = vIsPerspective == 0.0 ? gl_FragCoord.z : log2( vFragDepth ) * logDepthBufFC * 0.5;
#endif`, Bh = `#if defined( USE_LOGDEPTHBUF ) && defined( USE_LOGDEPTHBUF_EXT )
	uniform float logDepthBufFC;
	varying float vFragDepth;
	varying float vIsPerspective;
#endif`, zh = `#ifdef USE_LOGDEPTHBUF
	#ifdef USE_LOGDEPTHBUF_EXT
		varying float vFragDepth;
		varying float vIsPerspective;
	#else
		uniform float logDepthBufFC;
	#endif
#endif`, Hh = `#ifdef USE_LOGDEPTHBUF
	#ifdef USE_LOGDEPTHBUF_EXT
		vFragDepth = 1.0 + gl_Position.w;
		vIsPerspective = float( isPerspectiveMatrix( projectionMatrix ) );
	#else
		if ( isPerspectiveMatrix( projectionMatrix ) ) {
			gl_Position.z = log2( max( EPSILON, gl_Position.w + 1.0 ) ) * logDepthBufFC - 1.0;
			gl_Position.z *= gl_Position.w;
		}
	#endif
#endif`, Gh = `#ifdef USE_MAP
	diffuseColor *= texture2D( map, vMapUv );
#endif`, Vh = `#ifdef USE_MAP
	uniform sampler2D map;
#endif`, kh = `#if defined( USE_MAP ) || defined( USE_ALPHAMAP )
	#if defined( USE_POINTS_UV )
		vec2 uv = vUv;
	#else
		vec2 uv = ( uvTransform * vec3( gl_PointCoord.x, 1.0 - gl_PointCoord.y, 1 ) ).xy;
	#endif
#endif
#ifdef USE_MAP
	diffuseColor *= texture2D( map, uv );
#endif
#ifdef USE_ALPHAMAP
	diffuseColor.a *= texture2D( alphaMap, uv ).g;
#endif`, Wh = `#if defined( USE_POINTS_UV )
	varying vec2 vUv;
#else
	#if defined( USE_MAP ) || defined( USE_ALPHAMAP )
		uniform mat3 uvTransform;
	#endif
#endif
#ifdef USE_MAP
	uniform sampler2D map;
#endif
#ifdef USE_ALPHAMAP
	uniform sampler2D alphaMap;
#endif`, Xh = `float metalnessFactor = metalness;
#ifdef USE_METALNESSMAP
	vec4 texelMetalness = texture2D( metalnessMap, vMetalnessMapUv );
	metalnessFactor *= texelMetalness.b;
#endif`, Yh = `#ifdef USE_METALNESSMAP
	uniform sampler2D metalnessMap;
#endif`, qh = `#if defined( USE_MORPHCOLORS ) && defined( MORPHTARGETS_TEXTURE )
	vColor *= morphTargetBaseInfluence;
	for ( int i = 0; i < MORPHTARGETS_COUNT; i ++ ) {
		#if defined( USE_COLOR_ALPHA )
			if ( morphTargetInfluences[ i ] != 0.0 ) vColor += getMorph( gl_VertexID, i, 2 ) * morphTargetInfluences[ i ];
		#elif defined( USE_COLOR )
			if ( morphTargetInfluences[ i ] != 0.0 ) vColor += getMorph( gl_VertexID, i, 2 ).rgb * morphTargetInfluences[ i ];
		#endif
	}
#endif`, jh = `#ifdef USE_MORPHNORMALS
	objectNormal *= morphTargetBaseInfluence;
	#ifdef MORPHTARGETS_TEXTURE
		for ( int i = 0; i < MORPHTARGETS_COUNT; i ++ ) {
			if ( morphTargetInfluences[ i ] != 0.0 ) objectNormal += getMorph( gl_VertexID, i, 1 ).xyz * morphTargetInfluences[ i ];
		}
	#else
		objectNormal += morphNormal0 * morphTargetInfluences[ 0 ];
		objectNormal += morphNormal1 * morphTargetInfluences[ 1 ];
		objectNormal += morphNormal2 * morphTargetInfluences[ 2 ];
		objectNormal += morphNormal3 * morphTargetInfluences[ 3 ];
	#endif
#endif`, Zh = `#ifdef USE_MORPHTARGETS
	uniform float morphTargetBaseInfluence;
	#ifdef MORPHTARGETS_TEXTURE
		uniform float morphTargetInfluences[ MORPHTARGETS_COUNT ];
		uniform sampler2DArray morphTargetsTexture;
		uniform ivec2 morphTargetsTextureSize;
		vec4 getMorph( const in int vertexIndex, const in int morphTargetIndex, const in int offset ) {
			int texelIndex = vertexIndex * MORPHTARGETS_TEXTURE_STRIDE + offset;
			int y = texelIndex / morphTargetsTextureSize.x;
			int x = texelIndex - y * morphTargetsTextureSize.x;
			ivec3 morphUV = ivec3( x, y, morphTargetIndex );
			return texelFetch( morphTargetsTexture, morphUV, 0 );
		}
	#else
		#ifndef USE_MORPHNORMALS
			uniform float morphTargetInfluences[ 8 ];
		#else
			uniform float morphTargetInfluences[ 4 ];
		#endif
	#endif
#endif`, Kh = `#ifdef USE_MORPHTARGETS
	transformed *= morphTargetBaseInfluence;
	#ifdef MORPHTARGETS_TEXTURE
		for ( int i = 0; i < MORPHTARGETS_COUNT; i ++ ) {
			if ( morphTargetInfluences[ i ] != 0.0 ) transformed += getMorph( gl_VertexID, i, 0 ).xyz * morphTargetInfluences[ i ];
		}
	#else
		transformed += morphTarget0 * morphTargetInfluences[ 0 ];
		transformed += morphTarget1 * morphTargetInfluences[ 1 ];
		transformed += morphTarget2 * morphTargetInfluences[ 2 ];
		transformed += morphTarget3 * morphTargetInfluences[ 3 ];
		#ifndef USE_MORPHNORMALS
			transformed += morphTarget4 * morphTargetInfluences[ 4 ];
			transformed += morphTarget5 * morphTargetInfluences[ 5 ];
			transformed += morphTarget6 * morphTargetInfluences[ 6 ];
			transformed += morphTarget7 * morphTargetInfluences[ 7 ];
		#endif
	#endif
#endif`, Jh = `float faceDirection = gl_FrontFacing ? 1.0 : - 1.0;
#ifdef FLAT_SHADED
	vec3 fdx = dFdx( vViewPosition );
	vec3 fdy = dFdy( vViewPosition );
	vec3 normal = normalize( cross( fdx, fdy ) );
#else
	vec3 normal = normalize( vNormal );
	#ifdef DOUBLE_SIDED
		normal *= faceDirection;
	#endif
#endif
#if defined( USE_NORMALMAP_TANGENTSPACE ) || defined( USE_CLEARCOAT_NORMALMAP ) || defined( USE_ANISOTROPY )
	#ifdef USE_TANGENT
		mat3 tbn = mat3( normalize( vTangent ), normalize( vBitangent ), normal );
	#else
		mat3 tbn = getTangentFrame( - vViewPosition, normal,
		#if defined( USE_NORMALMAP )
			vNormalMapUv
		#elif defined( USE_CLEARCOAT_NORMALMAP )
			vClearcoatNormalMapUv
		#else
			vUv
		#endif
		);
	#endif
	#if defined( DOUBLE_SIDED ) && ! defined( FLAT_SHADED )
		tbn[0] *= faceDirection;
		tbn[1] *= faceDirection;
	#endif
#endif
#ifdef USE_CLEARCOAT_NORMALMAP
	#ifdef USE_TANGENT
		mat3 tbn2 = mat3( normalize( vTangent ), normalize( vBitangent ), normal );
	#else
		mat3 tbn2 = getTangentFrame( - vViewPosition, normal, vClearcoatNormalMapUv );
	#endif
	#if defined( DOUBLE_SIDED ) && ! defined( FLAT_SHADED )
		tbn2[0] *= faceDirection;
		tbn2[1] *= faceDirection;
	#endif
#endif
vec3 geometryNormal = normal;`, $h = `#ifdef USE_NORMALMAP_OBJECTSPACE
	normal = texture2D( normalMap, vNormalMapUv ).xyz * 2.0 - 1.0;
	#ifdef FLIP_SIDED
		normal = - normal;
	#endif
	#ifdef DOUBLE_SIDED
		normal = normal * faceDirection;
	#endif
	normal = normalize( normalMatrix * normal );
#elif defined( USE_NORMALMAP_TANGENTSPACE )
	vec3 mapN = texture2D( normalMap, vNormalMapUv ).xyz * 2.0 - 1.0;
	mapN.xy *= normalScale;
	normal = normalize( tbn * mapN );
#elif defined( USE_BUMPMAP )
	normal = perturbNormalArb( - vViewPosition, normal, dHdxy_fwd(), faceDirection );
#endif`, Qh = `#ifndef FLAT_SHADED
	varying vec3 vNormal;
	#ifdef USE_TANGENT
		varying vec3 vTangent;
		varying vec3 vBitangent;
	#endif
#endif`, eu = `#ifndef FLAT_SHADED
	varying vec3 vNormal;
	#ifdef USE_TANGENT
		varying vec3 vTangent;
		varying vec3 vBitangent;
	#endif
#endif`, tu = `#ifndef FLAT_SHADED
	vNormal = normalize( transformedNormal );
	#ifdef USE_TANGENT
		vTangent = normalize( transformedTangent );
		vBitangent = normalize( cross( vNormal, vTangent ) * tangent.w );
	#endif
#endif`, nu = `#ifdef USE_NORMALMAP
	uniform sampler2D normalMap;
	uniform vec2 normalScale;
#endif
#ifdef USE_NORMALMAP_OBJECTSPACE
	uniform mat3 normalMatrix;
#endif
#if ! defined ( USE_TANGENT ) && ( defined ( USE_NORMALMAP_TANGENTSPACE ) || defined ( USE_CLEARCOAT_NORMALMAP ) || defined( USE_ANISOTROPY ) )
	mat3 getTangentFrame( vec3 eye_pos, vec3 surf_norm, vec2 uv ) {
		vec3 q0 = dFdx( eye_pos.xyz );
		vec3 q1 = dFdy( eye_pos.xyz );
		vec2 st0 = dFdx( uv.st );
		vec2 st1 = dFdy( uv.st );
		vec3 N = surf_norm;
		vec3 q1perp = cross( q1, N );
		vec3 q0perp = cross( N, q0 );
		vec3 T = q1perp * st0.x + q0perp * st1.x;
		vec3 B = q1perp * st0.y + q0perp * st1.y;
		float det = max( dot( T, T ), dot( B, B ) );
		float scale = ( det == 0.0 ) ? 0.0 : inversesqrt( det );
		return mat3( T * scale, B * scale, N );
	}
#endif`, iu = `#ifdef USE_CLEARCOAT
	vec3 clearcoatNormal = geometryNormal;
#endif`, ru = `#ifdef USE_CLEARCOAT_NORMALMAP
	vec3 clearcoatMapN = texture2D( clearcoatNormalMap, vClearcoatNormalMapUv ).xyz * 2.0 - 1.0;
	clearcoatMapN.xy *= clearcoatNormalScale;
	clearcoatNormal = normalize( tbn2 * clearcoatMapN );
#endif`, su = `#ifdef USE_CLEARCOATMAP
	uniform sampler2D clearcoatMap;
#endif
#ifdef USE_CLEARCOAT_NORMALMAP
	uniform sampler2D clearcoatNormalMap;
	uniform vec2 clearcoatNormalScale;
#endif
#ifdef USE_CLEARCOAT_ROUGHNESSMAP
	uniform sampler2D clearcoatRoughnessMap;
#endif`, au = `#ifdef USE_IRIDESCENCEMAP
	uniform sampler2D iridescenceMap;
#endif
#ifdef USE_IRIDESCENCE_THICKNESSMAP
	uniform sampler2D iridescenceThicknessMap;
#endif`, ou = `#ifdef OPAQUE
diffuseColor.a = 1.0;
#endif
#ifdef USE_TRANSMISSION
diffuseColor.a *= material.transmissionAlpha;
#endif
gl_FragColor = vec4( outgoingLight, diffuseColor.a );`, lu = `vec3 packNormalToRGB( const in vec3 normal ) {
	return normalize( normal ) * 0.5 + 0.5;
}
vec3 unpackRGBToNormal( const in vec3 rgb ) {
	return 2.0 * rgb.xyz - 1.0;
}
const float PackUpscale = 256. / 255.;const float UnpackDownscale = 255. / 256.;
const vec3 PackFactors = vec3( 256. * 256. * 256., 256. * 256., 256. );
const vec4 UnpackFactors = UnpackDownscale / vec4( PackFactors, 1. );
const float ShiftRight8 = 1. / 256.;
vec4 packDepthToRGBA( const in float v ) {
	vec4 r = vec4( fract( v * PackFactors ), v );
	r.yzw -= r.xyz * ShiftRight8;	return r * PackUpscale;
}
float unpackRGBAToDepth( const in vec4 v ) {
	return dot( v, UnpackFactors );
}
vec2 packDepthToRG( in highp float v ) {
	return packDepthToRGBA( v ).yx;
}
float unpackRGToDepth( const in highp vec2 v ) {
	return unpackRGBAToDepth( vec4( v.xy, 0.0, 0.0 ) );
}
vec4 pack2HalfToRGBA( vec2 v ) {
	vec4 r = vec4( v.x, fract( v.x * 255.0 ), v.y, fract( v.y * 255.0 ) );
	return vec4( r.x - r.y / 255.0, r.y, r.z - r.w / 255.0, r.w );
}
vec2 unpackRGBATo2Half( vec4 v ) {
	return vec2( v.x + ( v.y / 255.0 ), v.z + ( v.w / 255.0 ) );
}
float viewZToOrthographicDepth( const in float viewZ, const in float near, const in float far ) {
	return ( viewZ + near ) / ( near - far );
}
float orthographicDepthToViewZ( const in float depth, const in float near, const in float far ) {
	return depth * ( near - far ) - near;
}
float viewZToPerspectiveDepth( const in float viewZ, const in float near, const in float far ) {
	return ( ( near + viewZ ) * far ) / ( ( far - near ) * viewZ );
}
float perspectiveDepthToViewZ( const in float depth, const in float near, const in float far ) {
	return ( near * far ) / ( ( far - near ) * depth - far );
}`, cu = `#ifdef PREMULTIPLIED_ALPHA
	gl_FragColor.rgb *= gl_FragColor.a;
#endif`, hu = `vec4 mvPosition = vec4( transformed, 1.0 );
#ifdef USE_INSTANCING
	mvPosition = instanceMatrix * mvPosition;
#endif
mvPosition = modelViewMatrix * mvPosition;
gl_Position = projectionMatrix * mvPosition;`, uu = `#ifdef DITHERING
	gl_FragColor.rgb = dithering( gl_FragColor.rgb );
#endif`, fu = `#ifdef DITHERING
	vec3 dithering( vec3 color ) {
		float grid_position = rand( gl_FragCoord.xy );
		vec3 dither_shift_RGB = vec3( 0.25 / 255.0, -0.25 / 255.0, 0.25 / 255.0 );
		dither_shift_RGB = mix( 2.0 * dither_shift_RGB, -2.0 * dither_shift_RGB, grid_position );
		return color + dither_shift_RGB;
	}
#endif`, du = `float roughnessFactor = roughness;
#ifdef USE_ROUGHNESSMAP
	vec4 texelRoughness = texture2D( roughnessMap, vRoughnessMapUv );
	roughnessFactor *= texelRoughness.g;
#endif`, pu = `#ifdef USE_ROUGHNESSMAP
	uniform sampler2D roughnessMap;
#endif`, mu = `#if NUM_SPOT_LIGHT_COORDS > 0
	varying vec4 vSpotLightCoord[ NUM_SPOT_LIGHT_COORDS ];
#endif
#if NUM_SPOT_LIGHT_MAPS > 0
	uniform sampler2D spotLightMap[ NUM_SPOT_LIGHT_MAPS ];
#endif
#ifdef USE_SHADOWMAP
	#if NUM_DIR_LIGHT_SHADOWS > 0
		uniform sampler2D directionalShadowMap[ NUM_DIR_LIGHT_SHADOWS ];
		varying vec4 vDirectionalShadowCoord[ NUM_DIR_LIGHT_SHADOWS ];
		struct DirectionalLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
		};
		uniform DirectionalLightShadow directionalLightShadows[ NUM_DIR_LIGHT_SHADOWS ];
	#endif
	#if NUM_SPOT_LIGHT_SHADOWS > 0
		uniform sampler2D spotShadowMap[ NUM_SPOT_LIGHT_SHADOWS ];
		struct SpotLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
		};
		uniform SpotLightShadow spotLightShadows[ NUM_SPOT_LIGHT_SHADOWS ];
	#endif
	#if NUM_POINT_LIGHT_SHADOWS > 0
		uniform sampler2D pointShadowMap[ NUM_POINT_LIGHT_SHADOWS ];
		varying vec4 vPointShadowCoord[ NUM_POINT_LIGHT_SHADOWS ];
		struct PointLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
			float shadowCameraNear;
			float shadowCameraFar;
		};
		uniform PointLightShadow pointLightShadows[ NUM_POINT_LIGHT_SHADOWS ];
	#endif
	float texture2DCompare( sampler2D depths, vec2 uv, float compare ) {
		return step( compare, unpackRGBAToDepth( texture2D( depths, uv ) ) );
	}
	vec2 texture2DDistribution( sampler2D shadow, vec2 uv ) {
		return unpackRGBATo2Half( texture2D( shadow, uv ) );
	}
	float VSMShadow (sampler2D shadow, vec2 uv, float compare ){
		float occlusion = 1.0;
		vec2 distribution = texture2DDistribution( shadow, uv );
		float hard_shadow = step( compare , distribution.x );
		if (hard_shadow != 1.0 ) {
			float distance = compare - distribution.x ;
			float variance = max( 0.00000, distribution.y * distribution.y );
			float softness_probability = variance / (variance + distance * distance );			softness_probability = clamp( ( softness_probability - 0.3 ) / ( 0.95 - 0.3 ), 0.0, 1.0 );			occlusion = clamp( max( hard_shadow, softness_probability ), 0.0, 1.0 );
		}
		return occlusion;
	}
	float getShadow( sampler2D shadowMap, vec2 shadowMapSize, float shadowBias, float shadowRadius, vec4 shadowCoord ) {
		float shadow = 1.0;
		shadowCoord.xyz /= shadowCoord.w;
		shadowCoord.z += shadowBias;
		bool inFrustum = shadowCoord.x >= 0.0 && shadowCoord.x <= 1.0 && shadowCoord.y >= 0.0 && shadowCoord.y <= 1.0;
		bool frustumTest = inFrustum && shadowCoord.z <= 1.0;
		if ( frustumTest ) {
		#if defined( SHADOWMAP_TYPE_PCF )
			vec2 texelSize = vec2( 1.0 ) / shadowMapSize;
			float dx0 = - texelSize.x * shadowRadius;
			float dy0 = - texelSize.y * shadowRadius;
			float dx1 = + texelSize.x * shadowRadius;
			float dy1 = + texelSize.y * shadowRadius;
			float dx2 = dx0 / 2.0;
			float dy2 = dy0 / 2.0;
			float dx3 = dx1 / 2.0;
			float dy3 = dy1 / 2.0;
			shadow = (
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, dy0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, dy0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, dy2 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy2 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, dy2 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, 0.0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, 0.0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy, shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, 0.0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, 0.0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx2, dy3 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy3 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx3, dy3 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx0, dy1 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( 0.0, dy1 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, shadowCoord.xy + vec2( dx1, dy1 ), shadowCoord.z )
			) * ( 1.0 / 17.0 );
		#elif defined( SHADOWMAP_TYPE_PCF_SOFT )
			vec2 texelSize = vec2( 1.0 ) / shadowMapSize;
			float dx = texelSize.x;
			float dy = texelSize.y;
			vec2 uv = shadowCoord.xy;
			vec2 f = fract( uv * shadowMapSize + 0.5 );
			uv -= f * texelSize;
			shadow = (
				texture2DCompare( shadowMap, uv, shadowCoord.z ) +
				texture2DCompare( shadowMap, uv + vec2( dx, 0.0 ), shadowCoord.z ) +
				texture2DCompare( shadowMap, uv + vec2( 0.0, dy ), shadowCoord.z ) +
				texture2DCompare( shadowMap, uv + texelSize, shadowCoord.z ) +
				mix( texture2DCompare( shadowMap, uv + vec2( -dx, 0.0 ), shadowCoord.z ),
					 texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, 0.0 ), shadowCoord.z ),
					 f.x ) +
				mix( texture2DCompare( shadowMap, uv + vec2( -dx, dy ), shadowCoord.z ),
					 texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, dy ), shadowCoord.z ),
					 f.x ) +
				mix( texture2DCompare( shadowMap, uv + vec2( 0.0, -dy ), shadowCoord.z ),
					 texture2DCompare( shadowMap, uv + vec2( 0.0, 2.0 * dy ), shadowCoord.z ),
					 f.y ) +
				mix( texture2DCompare( shadowMap, uv + vec2( dx, -dy ), shadowCoord.z ),
					 texture2DCompare( shadowMap, uv + vec2( dx, 2.0 * dy ), shadowCoord.z ),
					 f.y ) +
				mix( mix( texture2DCompare( shadowMap, uv + vec2( -dx, -dy ), shadowCoord.z ),
						  texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, -dy ), shadowCoord.z ),
						  f.x ),
					 mix( texture2DCompare( shadowMap, uv + vec2( -dx, 2.0 * dy ), shadowCoord.z ),
						  texture2DCompare( shadowMap, uv + vec2( 2.0 * dx, 2.0 * dy ), shadowCoord.z ),
						  f.x ),
					 f.y )
			) * ( 1.0 / 9.0 );
		#elif defined( SHADOWMAP_TYPE_VSM )
			shadow = VSMShadow( shadowMap, shadowCoord.xy, shadowCoord.z );
		#else
			shadow = texture2DCompare( shadowMap, shadowCoord.xy, shadowCoord.z );
		#endif
		}
		return shadow;
	}
	vec2 cubeToUV( vec3 v, float texelSizeY ) {
		vec3 absV = abs( v );
		float scaleToCube = 1.0 / max( absV.x, max( absV.y, absV.z ) );
		absV *= scaleToCube;
		v *= scaleToCube * ( 1.0 - 2.0 * texelSizeY );
		vec2 planar = v.xy;
		float almostATexel = 1.5 * texelSizeY;
		float almostOne = 1.0 - almostATexel;
		if ( absV.z >= almostOne ) {
			if ( v.z > 0.0 )
				planar.x = 4.0 - v.x;
		} else if ( absV.x >= almostOne ) {
			float signX = sign( v.x );
			planar.x = v.z * signX + 2.0 * signX;
		} else if ( absV.y >= almostOne ) {
			float signY = sign( v.y );
			planar.x = v.x + 2.0 * signY + 2.0;
			planar.y = v.z * signY - 2.0;
		}
		return vec2( 0.125, 0.25 ) * planar + vec2( 0.375, 0.75 );
	}
	float getPointShadow( sampler2D shadowMap, vec2 shadowMapSize, float shadowBias, float shadowRadius, vec4 shadowCoord, float shadowCameraNear, float shadowCameraFar ) {
		vec2 texelSize = vec2( 1.0 ) / ( shadowMapSize * vec2( 4.0, 2.0 ) );
		vec3 lightToPosition = shadowCoord.xyz;
		float dp = ( length( lightToPosition ) - shadowCameraNear ) / ( shadowCameraFar - shadowCameraNear );		dp += shadowBias;
		vec3 bd3D = normalize( lightToPosition );
		#if defined( SHADOWMAP_TYPE_PCF ) || defined( SHADOWMAP_TYPE_PCF_SOFT ) || defined( SHADOWMAP_TYPE_VSM )
			vec2 offset = vec2( - 1, 1 ) * shadowRadius * texelSize.y;
			return (
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.xyy, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.yyy, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.xyx, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.yyx, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.xxy, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.yxy, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.xxx, texelSize.y ), dp ) +
				texture2DCompare( shadowMap, cubeToUV( bd3D + offset.yxx, texelSize.y ), dp )
			) * ( 1.0 / 9.0 );
		#else
			return texture2DCompare( shadowMap, cubeToUV( bd3D, texelSize.y ), dp );
		#endif
	}
#endif`, gu = `#if NUM_SPOT_LIGHT_COORDS > 0
	uniform mat4 spotLightMatrix[ NUM_SPOT_LIGHT_COORDS ];
	varying vec4 vSpotLightCoord[ NUM_SPOT_LIGHT_COORDS ];
#endif
#ifdef USE_SHADOWMAP
	#if NUM_DIR_LIGHT_SHADOWS > 0
		uniform mat4 directionalShadowMatrix[ NUM_DIR_LIGHT_SHADOWS ];
		varying vec4 vDirectionalShadowCoord[ NUM_DIR_LIGHT_SHADOWS ];
		struct DirectionalLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
		};
		uniform DirectionalLightShadow directionalLightShadows[ NUM_DIR_LIGHT_SHADOWS ];
	#endif
	#if NUM_SPOT_LIGHT_SHADOWS > 0
		struct SpotLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
		};
		uniform SpotLightShadow spotLightShadows[ NUM_SPOT_LIGHT_SHADOWS ];
	#endif
	#if NUM_POINT_LIGHT_SHADOWS > 0
		uniform mat4 pointShadowMatrix[ NUM_POINT_LIGHT_SHADOWS ];
		varying vec4 vPointShadowCoord[ NUM_POINT_LIGHT_SHADOWS ];
		struct PointLightShadow {
			float shadowBias;
			float shadowNormalBias;
			float shadowRadius;
			vec2 shadowMapSize;
			float shadowCameraNear;
			float shadowCameraFar;
		};
		uniform PointLightShadow pointLightShadows[ NUM_POINT_LIGHT_SHADOWS ];
	#endif
#endif`, _u = `#if ( defined( USE_SHADOWMAP ) && ( NUM_DIR_LIGHT_SHADOWS > 0 || NUM_POINT_LIGHT_SHADOWS > 0 ) ) || ( NUM_SPOT_LIGHT_COORDS > 0 )
	vec3 shadowWorldNormal = inverseTransformDirection( transformedNormal, viewMatrix );
	vec4 shadowWorldPosition;
#endif
#if defined( USE_SHADOWMAP )
	#if NUM_DIR_LIGHT_SHADOWS > 0
		#pragma unroll_loop_start
		for ( int i = 0; i < NUM_DIR_LIGHT_SHADOWS; i ++ ) {
			shadowWorldPosition = worldPosition + vec4( shadowWorldNormal * directionalLightShadows[ i ].shadowNormalBias, 0 );
			vDirectionalShadowCoord[ i ] = directionalShadowMatrix[ i ] * shadowWorldPosition;
		}
		#pragma unroll_loop_end
	#endif
	#if NUM_POINT_LIGHT_SHADOWS > 0
		#pragma unroll_loop_start
		for ( int i = 0; i < NUM_POINT_LIGHT_SHADOWS; i ++ ) {
			shadowWorldPosition = worldPosition + vec4( shadowWorldNormal * pointLightShadows[ i ].shadowNormalBias, 0 );
			vPointShadowCoord[ i ] = pointShadowMatrix[ i ] * shadowWorldPosition;
		}
		#pragma unroll_loop_end
	#endif
#endif
#if NUM_SPOT_LIGHT_COORDS > 0
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_SPOT_LIGHT_COORDS; i ++ ) {
		shadowWorldPosition = worldPosition;
		#if ( defined( USE_SHADOWMAP ) && UNROLLED_LOOP_INDEX < NUM_SPOT_LIGHT_SHADOWS )
			shadowWorldPosition.xyz += shadowWorldNormal * spotLightShadows[ i ].shadowNormalBias;
		#endif
		vSpotLightCoord[ i ] = spotLightMatrix[ i ] * shadowWorldPosition;
	}
	#pragma unroll_loop_end
#endif`, vu = `float getShadowMask() {
	float shadow = 1.0;
	#ifdef USE_SHADOWMAP
	#if NUM_DIR_LIGHT_SHADOWS > 0
	DirectionalLightShadow directionalLight;
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_DIR_LIGHT_SHADOWS; i ++ ) {
		directionalLight = directionalLightShadows[ i ];
		shadow *= receiveShadow ? getShadow( directionalShadowMap[ i ], directionalLight.shadowMapSize, directionalLight.shadowBias, directionalLight.shadowRadius, vDirectionalShadowCoord[ i ] ) : 1.0;
	}
	#pragma unroll_loop_end
	#endif
	#if NUM_SPOT_LIGHT_SHADOWS > 0
	SpotLightShadow spotLight;
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_SPOT_LIGHT_SHADOWS; i ++ ) {
		spotLight = spotLightShadows[ i ];
		shadow *= receiveShadow ? getShadow( spotShadowMap[ i ], spotLight.shadowMapSize, spotLight.shadowBias, spotLight.shadowRadius, vSpotLightCoord[ i ] ) : 1.0;
	}
	#pragma unroll_loop_end
	#endif
	#if NUM_POINT_LIGHT_SHADOWS > 0
	PointLightShadow pointLight;
	#pragma unroll_loop_start
	for ( int i = 0; i < NUM_POINT_LIGHT_SHADOWS; i ++ ) {
		pointLight = pointLightShadows[ i ];
		shadow *= receiveShadow ? getPointShadow( pointShadowMap[ i ], pointLight.shadowMapSize, pointLight.shadowBias, pointLight.shadowRadius, vPointShadowCoord[ i ], pointLight.shadowCameraNear, pointLight.shadowCameraFar ) : 1.0;
	}
	#pragma unroll_loop_end
	#endif
	#endif
	return shadow;
}`, xu = `#ifdef USE_SKINNING
	mat4 boneMatX = getBoneMatrix( skinIndex.x );
	mat4 boneMatY = getBoneMatrix( skinIndex.y );
	mat4 boneMatZ = getBoneMatrix( skinIndex.z );
	mat4 boneMatW = getBoneMatrix( skinIndex.w );
#endif`, Mu = `#ifdef USE_SKINNING
	uniform mat4 bindMatrix;
	uniform mat4 bindMatrixInverse;
	uniform highp sampler2D boneTexture;
	uniform int boneTextureSize;
	mat4 getBoneMatrix( const in float i ) {
		float j = i * 4.0;
		float x = mod( j, float( boneTextureSize ) );
		float y = floor( j / float( boneTextureSize ) );
		float dx = 1.0 / float( boneTextureSize );
		float dy = 1.0 / float( boneTextureSize );
		y = dy * ( y + 0.5 );
		vec4 v1 = texture2D( boneTexture, vec2( dx * ( x + 0.5 ), y ) );
		vec4 v2 = texture2D( boneTexture, vec2( dx * ( x + 1.5 ), y ) );
		vec4 v3 = texture2D( boneTexture, vec2( dx * ( x + 2.5 ), y ) );
		vec4 v4 = texture2D( boneTexture, vec2( dx * ( x + 3.5 ), y ) );
		mat4 bone = mat4( v1, v2, v3, v4 );
		return bone;
	}
#endif`, Su = `#ifdef USE_SKINNING
	vec4 skinVertex = bindMatrix * vec4( transformed, 1.0 );
	vec4 skinned = vec4( 0.0 );
	skinned += boneMatX * skinVertex * skinWeight.x;
	skinned += boneMatY * skinVertex * skinWeight.y;
	skinned += boneMatZ * skinVertex * skinWeight.z;
	skinned += boneMatW * skinVertex * skinWeight.w;
	transformed = ( bindMatrixInverse * skinned ).xyz;
#endif`, Eu = `#ifdef USE_SKINNING
	mat4 skinMatrix = mat4( 0.0 );
	skinMatrix += skinWeight.x * boneMatX;
	skinMatrix += skinWeight.y * boneMatY;
	skinMatrix += skinWeight.z * boneMatZ;
	skinMatrix += skinWeight.w * boneMatW;
	skinMatrix = bindMatrixInverse * skinMatrix * bindMatrix;
	objectNormal = vec4( skinMatrix * vec4( objectNormal, 0.0 ) ).xyz;
	#ifdef USE_TANGENT
		objectTangent = vec4( skinMatrix * vec4( objectTangent, 0.0 ) ).xyz;
	#endif
#endif`, yu = `float specularStrength;
#ifdef USE_SPECULARMAP
	vec4 texelSpecular = texture2D( specularMap, vSpecularMapUv );
	specularStrength = texelSpecular.r;
#else
	specularStrength = 1.0;
#endif`, Tu = `#ifdef USE_SPECULARMAP
	uniform sampler2D specularMap;
#endif`, Au = `#if defined( TONE_MAPPING )
	gl_FragColor.rgb = toneMapping( gl_FragColor.rgb );
#endif`, bu = `#ifndef saturate
#define saturate( a ) clamp( a, 0.0, 1.0 )
#endif
uniform float toneMappingExposure;
vec3 LinearToneMapping( vec3 color ) {
	return saturate( toneMappingExposure * color );
}
vec3 ReinhardToneMapping( vec3 color ) {
	color *= toneMappingExposure;
	return saturate( color / ( vec3( 1.0 ) + color ) );
}
vec3 OptimizedCineonToneMapping( vec3 color ) {
	color *= toneMappingExposure;
	color = max( vec3( 0.0 ), color - 0.004 );
	return pow( ( color * ( 6.2 * color + 0.5 ) ) / ( color * ( 6.2 * color + 1.7 ) + 0.06 ), vec3( 2.2 ) );
}
vec3 RRTAndODTFit( vec3 v ) {
	vec3 a = v * ( v + 0.0245786 ) - 0.000090537;
	vec3 b = v * ( 0.983729 * v + 0.4329510 ) + 0.238081;
	return a / b;
}
vec3 ACESFilmicToneMapping( vec3 color ) {
	const mat3 ACESInputMat = mat3(
		vec3( 0.59719, 0.07600, 0.02840 ),		vec3( 0.35458, 0.90834, 0.13383 ),
		vec3( 0.04823, 0.01566, 0.83777 )
	);
	const mat3 ACESOutputMat = mat3(
		vec3(  1.60475, -0.10208, -0.00327 ),		vec3( -0.53108,  1.10813, -0.07276 ),
		vec3( -0.07367, -0.00605,  1.07602 )
	);
	color *= toneMappingExposure / 0.6;
	color = ACESInputMat * color;
	color = RRTAndODTFit( color );
	color = ACESOutputMat * color;
	return saturate( color );
}
vec3 CustomToneMapping( vec3 color ) { return color; }`, wu = `#ifdef USE_TRANSMISSION
	material.transmission = transmission;
	material.transmissionAlpha = 1.0;
	material.thickness = thickness;
	material.attenuationDistance = attenuationDistance;
	material.attenuationColor = attenuationColor;
	#ifdef USE_TRANSMISSIONMAP
		material.transmission *= texture2D( transmissionMap, vTransmissionMapUv ).r;
	#endif
	#ifdef USE_THICKNESSMAP
		material.thickness *= texture2D( thicknessMap, vThicknessMapUv ).g;
	#endif
	vec3 pos = vWorldPosition;
	vec3 v = normalize( cameraPosition - pos );
	vec3 n = inverseTransformDirection( normal, viewMatrix );
	vec4 transmitted = getIBLVolumeRefraction(
		n, v, material.roughness, material.diffuseColor, material.specularColor, material.specularF90,
		pos, modelMatrix, viewMatrix, projectionMatrix, material.ior, material.thickness,
		material.attenuationColor, material.attenuationDistance );
	material.transmissionAlpha = mix( material.transmissionAlpha, transmitted.a, material.transmission );
	totalDiffuse = mix( totalDiffuse, transmitted.rgb, material.transmission );
#endif`, Ru = `#ifdef USE_TRANSMISSION
	uniform float transmission;
	uniform float thickness;
	uniform float attenuationDistance;
	uniform vec3 attenuationColor;
	#ifdef USE_TRANSMISSIONMAP
		uniform sampler2D transmissionMap;
	#endif
	#ifdef USE_THICKNESSMAP
		uniform sampler2D thicknessMap;
	#endif
	uniform vec2 transmissionSamplerSize;
	uniform sampler2D transmissionSamplerMap;
	uniform mat4 modelMatrix;
	uniform mat4 projectionMatrix;
	varying vec3 vWorldPosition;
	float w0( float a ) {
		return ( 1.0 / 6.0 ) * ( a * ( a * ( - a + 3.0 ) - 3.0 ) + 1.0 );
	}
	float w1( float a ) {
		return ( 1.0 / 6.0 ) * ( a *  a * ( 3.0 * a - 6.0 ) + 4.0 );
	}
	float w2( float a ){
		return ( 1.0 / 6.0 ) * ( a * ( a * ( - 3.0 * a + 3.0 ) + 3.0 ) + 1.0 );
	}
	float w3( float a ) {
		return ( 1.0 / 6.0 ) * ( a * a * a );
	}
	float g0( float a ) {
		return w0( a ) + w1( a );
	}
	float g1( float a ) {
		return w2( a ) + w3( a );
	}
	float h0( float a ) {
		return - 1.0 + w1( a ) / ( w0( a ) + w1( a ) );
	}
	float h1( float a ) {
		return 1.0 + w3( a ) / ( w2( a ) + w3( a ) );
	}
	vec4 bicubic( sampler2D tex, vec2 uv, vec4 texelSize, float lod ) {
		uv = uv * texelSize.zw + 0.5;
		vec2 iuv = floor( uv );
		vec2 fuv = fract( uv );
		float g0x = g0( fuv.x );
		float g1x = g1( fuv.x );
		float h0x = h0( fuv.x );
		float h1x = h1( fuv.x );
		float h0y = h0( fuv.y );
		float h1y = h1( fuv.y );
		vec2 p0 = ( vec2( iuv.x + h0x, iuv.y + h0y ) - 0.5 ) * texelSize.xy;
		vec2 p1 = ( vec2( iuv.x + h1x, iuv.y + h0y ) - 0.5 ) * texelSize.xy;
		vec2 p2 = ( vec2( iuv.x + h0x, iuv.y + h1y ) - 0.5 ) * texelSize.xy;
		vec2 p3 = ( vec2( iuv.x + h1x, iuv.y + h1y ) - 0.5 ) * texelSize.xy;
		return g0( fuv.y ) * ( g0x * textureLod( tex, p0, lod ) + g1x * textureLod( tex, p1, lod ) ) +
			g1( fuv.y ) * ( g0x * textureLod( tex, p2, lod ) + g1x * textureLod( tex, p3, lod ) );
	}
	vec4 textureBicubic( sampler2D sampler, vec2 uv, float lod ) {
		vec2 fLodSize = vec2( textureSize( sampler, int( lod ) ) );
		vec2 cLodSize = vec2( textureSize( sampler, int( lod + 1.0 ) ) );
		vec2 fLodSizeInv = 1.0 / fLodSize;
		vec2 cLodSizeInv = 1.0 / cLodSize;
		vec4 fSample = bicubic( sampler, uv, vec4( fLodSizeInv, fLodSize ), floor( lod ) );
		vec4 cSample = bicubic( sampler, uv, vec4( cLodSizeInv, cLodSize ), ceil( lod ) );
		return mix( fSample, cSample, fract( lod ) );
	}
	vec3 getVolumeTransmissionRay( const in vec3 n, const in vec3 v, const in float thickness, const in float ior, const in mat4 modelMatrix ) {
		vec3 refractionVector = refract( - v, normalize( n ), 1.0 / ior );
		vec3 modelScale;
		modelScale.x = length( vec3( modelMatrix[ 0 ].xyz ) );
		modelScale.y = length( vec3( modelMatrix[ 1 ].xyz ) );
		modelScale.z = length( vec3( modelMatrix[ 2 ].xyz ) );
		return normalize( refractionVector ) * thickness * modelScale;
	}
	float applyIorToRoughness( const in float roughness, const in float ior ) {
		return roughness * clamp( ior * 2.0 - 2.0, 0.0, 1.0 );
	}
	vec4 getTransmissionSample( const in vec2 fragCoord, const in float roughness, const in float ior ) {
		float lod = log2( transmissionSamplerSize.x ) * applyIorToRoughness( roughness, ior );
		return textureBicubic( transmissionSamplerMap, fragCoord.xy, lod );
	}
	vec3 volumeAttenuation( const in float transmissionDistance, const in vec3 attenuationColor, const in float attenuationDistance ) {
		if ( isinf( attenuationDistance ) ) {
			return vec3( 1.0 );
		} else {
			vec3 attenuationCoefficient = -log( attenuationColor ) / attenuationDistance;
			vec3 transmittance = exp( - attenuationCoefficient * transmissionDistance );			return transmittance;
		}
	}
	vec4 getIBLVolumeRefraction( const in vec3 n, const in vec3 v, const in float roughness, const in vec3 diffuseColor,
		const in vec3 specularColor, const in float specularF90, const in vec3 position, const in mat4 modelMatrix,
		const in mat4 viewMatrix, const in mat4 projMatrix, const in float ior, const in float thickness,
		const in vec3 attenuationColor, const in float attenuationDistance ) {
		vec3 transmissionRay = getVolumeTransmissionRay( n, v, thickness, ior, modelMatrix );
		vec3 refractedRayExit = position + transmissionRay;
		vec4 ndcPos = projMatrix * viewMatrix * vec4( refractedRayExit, 1.0 );
		vec2 refractionCoords = ndcPos.xy / ndcPos.w;
		refractionCoords += 1.0;
		refractionCoords /= 2.0;
		vec4 transmittedLight = getTransmissionSample( refractionCoords, roughness, ior );
		vec3 transmittance = diffuseColor * volumeAttenuation( length( transmissionRay ), attenuationColor, attenuationDistance );
		vec3 attenuatedColor = transmittance * transmittedLight.rgb;
		vec3 F = EnvironmentBRDF( n, v, specularColor, specularF90, roughness );
		float transmittanceFactor = ( transmittance.r + transmittance.g + transmittance.b ) / 3.0;
		return vec4( ( 1.0 - F ) * attenuatedColor, 1.0 - ( 1.0 - transmittedLight.a ) * transmittanceFactor );
	}
#endif`, Cu = `#if defined( USE_UV ) || defined( USE_ANISOTROPY )
	varying vec2 vUv;
#endif
#ifdef USE_MAP
	varying vec2 vMapUv;
#endif
#ifdef USE_ALPHAMAP
	varying vec2 vAlphaMapUv;
#endif
#ifdef USE_LIGHTMAP
	varying vec2 vLightMapUv;
#endif
#ifdef USE_AOMAP
	varying vec2 vAoMapUv;
#endif
#ifdef USE_BUMPMAP
	varying vec2 vBumpMapUv;
#endif
#ifdef USE_NORMALMAP
	varying vec2 vNormalMapUv;
#endif
#ifdef USE_EMISSIVEMAP
	varying vec2 vEmissiveMapUv;
#endif
#ifdef USE_METALNESSMAP
	varying vec2 vMetalnessMapUv;
#endif
#ifdef USE_ROUGHNESSMAP
	varying vec2 vRoughnessMapUv;
#endif
#ifdef USE_ANISOTROPYMAP
	varying vec2 vAnisotropyMapUv;
#endif
#ifdef USE_CLEARCOATMAP
	varying vec2 vClearcoatMapUv;
#endif
#ifdef USE_CLEARCOAT_NORMALMAP
	varying vec2 vClearcoatNormalMapUv;
#endif
#ifdef USE_CLEARCOAT_ROUGHNESSMAP
	varying vec2 vClearcoatRoughnessMapUv;
#endif
#ifdef USE_IRIDESCENCEMAP
	varying vec2 vIridescenceMapUv;
#endif
#ifdef USE_IRIDESCENCE_THICKNESSMAP
	varying vec2 vIridescenceThicknessMapUv;
#endif
#ifdef USE_SHEEN_COLORMAP
	varying vec2 vSheenColorMapUv;
#endif
#ifdef USE_SHEEN_ROUGHNESSMAP
	varying vec2 vSheenRoughnessMapUv;
#endif
#ifdef USE_SPECULARMAP
	varying vec2 vSpecularMapUv;
#endif
#ifdef USE_SPECULAR_COLORMAP
	varying vec2 vSpecularColorMapUv;
#endif
#ifdef USE_SPECULAR_INTENSITYMAP
	varying vec2 vSpecularIntensityMapUv;
#endif
#ifdef USE_TRANSMISSIONMAP
	uniform mat3 transmissionMapTransform;
	varying vec2 vTransmissionMapUv;
#endif
#ifdef USE_THICKNESSMAP
	uniform mat3 thicknessMapTransform;
	varying vec2 vThicknessMapUv;
#endif`, Pu = `#if defined( USE_UV ) || defined( USE_ANISOTROPY )
	varying vec2 vUv;
#endif
#ifdef USE_MAP
	uniform mat3 mapTransform;
	varying vec2 vMapUv;
#endif
#ifdef USE_ALPHAMAP
	uniform mat3 alphaMapTransform;
	varying vec2 vAlphaMapUv;
#endif
#ifdef USE_LIGHTMAP
	uniform mat3 lightMapTransform;
	varying vec2 vLightMapUv;
#endif
#ifdef USE_AOMAP
	uniform mat3 aoMapTransform;
	varying vec2 vAoMapUv;
#endif
#ifdef USE_BUMPMAP
	uniform mat3 bumpMapTransform;
	varying vec2 vBumpMapUv;
#endif
#ifdef USE_NORMALMAP
	uniform mat3 normalMapTransform;
	varying vec2 vNormalMapUv;
#endif
#ifdef USE_DISPLACEMENTMAP
	uniform mat3 displacementMapTransform;
	varying vec2 vDisplacementMapUv;
#endif
#ifdef USE_EMISSIVEMAP
	uniform mat3 emissiveMapTransform;
	varying vec2 vEmissiveMapUv;
#endif
#ifdef USE_METALNESSMAP
	uniform mat3 metalnessMapTransform;
	varying vec2 vMetalnessMapUv;
#endif
#ifdef USE_ROUGHNESSMAP
	uniform mat3 roughnessMapTransform;
	varying vec2 vRoughnessMapUv;
#endif
#ifdef USE_ANISOTROPYMAP
	uniform mat3 anisotropyMapTransform;
	varying vec2 vAnisotropyMapUv;
#endif
#ifdef USE_CLEARCOATMAP
	uniform mat3 clearcoatMapTransform;
	varying vec2 vClearcoatMapUv;
#endif
#ifdef USE_CLEARCOAT_NORMALMAP
	uniform mat3 clearcoatNormalMapTransform;
	varying vec2 vClearcoatNormalMapUv;
#endif
#ifdef USE_CLEARCOAT_ROUGHNESSMAP
	uniform mat3 clearcoatRoughnessMapTransform;
	varying vec2 vClearcoatRoughnessMapUv;
#endif
#ifdef USE_SHEEN_COLORMAP
	uniform mat3 sheenColorMapTransform;
	varying vec2 vSheenColorMapUv;
#endif
#ifdef USE_SHEEN_ROUGHNESSMAP
	uniform mat3 sheenRoughnessMapTransform;
	varying vec2 vSheenRoughnessMapUv;
#endif
#ifdef USE_IRIDESCENCEMAP
	uniform mat3 iridescenceMapTransform;
	varying vec2 vIridescenceMapUv;
#endif
#ifdef USE_IRIDESCENCE_THICKNESSMAP
	uniform mat3 iridescenceThicknessMapTransform;
	varying vec2 vIridescenceThicknessMapUv;
#endif
#ifdef USE_SPECULARMAP
	uniform mat3 specularMapTransform;
	varying vec2 vSpecularMapUv;
#endif
#ifdef USE_SPECULAR_COLORMAP
	uniform mat3 specularColorMapTransform;
	varying vec2 vSpecularColorMapUv;
#endif
#ifdef USE_SPECULAR_INTENSITYMAP
	uniform mat3 specularIntensityMapTransform;
	varying vec2 vSpecularIntensityMapUv;
#endif
#ifdef USE_TRANSMISSIONMAP
	uniform mat3 transmissionMapTransform;
	varying vec2 vTransmissionMapUv;
#endif
#ifdef USE_THICKNESSMAP
	uniform mat3 thicknessMapTransform;
	varying vec2 vThicknessMapUv;
#endif`, Lu = `#if defined( USE_UV ) || defined( USE_ANISOTROPY )
	vUv = vec3( uv, 1 ).xy;
#endif
#ifdef USE_MAP
	vMapUv = ( mapTransform * vec3( MAP_UV, 1 ) ).xy;
#endif
#ifdef USE_ALPHAMAP
	vAlphaMapUv = ( alphaMapTransform * vec3( ALPHAMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_LIGHTMAP
	vLightMapUv = ( lightMapTransform * vec3( LIGHTMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_AOMAP
	vAoMapUv = ( aoMapTransform * vec3( AOMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_BUMPMAP
	vBumpMapUv = ( bumpMapTransform * vec3( BUMPMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_NORMALMAP
	vNormalMapUv = ( normalMapTransform * vec3( NORMALMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_DISPLACEMENTMAP
	vDisplacementMapUv = ( displacementMapTransform * vec3( DISPLACEMENTMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_EMISSIVEMAP
	vEmissiveMapUv = ( emissiveMapTransform * vec3( EMISSIVEMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_METALNESSMAP
	vMetalnessMapUv = ( metalnessMapTransform * vec3( METALNESSMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_ROUGHNESSMAP
	vRoughnessMapUv = ( roughnessMapTransform * vec3( ROUGHNESSMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_ANISOTROPYMAP
	vAnisotropyMapUv = ( anisotropyMapTransform * vec3( ANISOTROPYMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_CLEARCOATMAP
	vClearcoatMapUv = ( clearcoatMapTransform * vec3( CLEARCOATMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_CLEARCOAT_NORMALMAP
	vClearcoatNormalMapUv = ( clearcoatNormalMapTransform * vec3( CLEARCOAT_NORMALMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_CLEARCOAT_ROUGHNESSMAP
	vClearcoatRoughnessMapUv = ( clearcoatRoughnessMapTransform * vec3( CLEARCOAT_ROUGHNESSMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_IRIDESCENCEMAP
	vIridescenceMapUv = ( iridescenceMapTransform * vec3( IRIDESCENCEMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_IRIDESCENCE_THICKNESSMAP
	vIridescenceThicknessMapUv = ( iridescenceThicknessMapTransform * vec3( IRIDESCENCE_THICKNESSMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_SHEEN_COLORMAP
	vSheenColorMapUv = ( sheenColorMapTransform * vec3( SHEEN_COLORMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_SHEEN_ROUGHNESSMAP
	vSheenRoughnessMapUv = ( sheenRoughnessMapTransform * vec3( SHEEN_ROUGHNESSMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_SPECULARMAP
	vSpecularMapUv = ( specularMapTransform * vec3( SPECULARMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_SPECULAR_COLORMAP
	vSpecularColorMapUv = ( specularColorMapTransform * vec3( SPECULAR_COLORMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_SPECULAR_INTENSITYMAP
	vSpecularIntensityMapUv = ( specularIntensityMapTransform * vec3( SPECULAR_INTENSITYMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_TRANSMISSIONMAP
	vTransmissionMapUv = ( transmissionMapTransform * vec3( TRANSMISSIONMAP_UV, 1 ) ).xy;
#endif
#ifdef USE_THICKNESSMAP
	vThicknessMapUv = ( thicknessMapTransform * vec3( THICKNESSMAP_UV, 1 ) ).xy;
#endif`, Uu = `#if defined( USE_ENVMAP ) || defined( DISTANCE ) || defined ( USE_SHADOWMAP ) || defined ( USE_TRANSMISSION ) || NUM_SPOT_LIGHT_COORDS > 0
	vec4 worldPosition = vec4( transformed, 1.0 );
	#ifdef USE_INSTANCING
		worldPosition = instanceMatrix * worldPosition;
	#endif
	worldPosition = modelMatrix * worldPosition;
#endif`;
const Du = `varying vec2 vUv;
uniform mat3 uvTransform;
void main() {
	vUv = ( uvTransform * vec3( uv, 1 ) ).xy;
	gl_Position = vec4( position.xy, 1.0, 1.0 );
}`, Iu = `uniform sampler2D t2D;
uniform float backgroundIntensity;
varying vec2 vUv;
void main() {
	vec4 texColor = texture2D( t2D, vUv );
	texColor.rgb *= backgroundIntensity;
	gl_FragColor = texColor;
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
}`, Nu = `varying vec3 vWorldDirection;
#include <common>
void main() {
	vWorldDirection = transformDirection( position, modelMatrix );
	#include <begin_vertex>
	#include <project_vertex>
	gl_Position.z = gl_Position.w;
}`, Ou = `#ifdef ENVMAP_TYPE_CUBE
	uniform samplerCube envMap;
#elif defined( ENVMAP_TYPE_CUBE_UV )
	uniform sampler2D envMap;
#endif
uniform float flipEnvMap;
uniform float backgroundBlurriness;
uniform float backgroundIntensity;
varying vec3 vWorldDirection;
#include <cube_uv_reflection_fragment>
void main() {
	#ifdef ENVMAP_TYPE_CUBE
		vec4 texColor = textureCube( envMap, vec3( flipEnvMap * vWorldDirection.x, vWorldDirection.yz ) );
	#elif defined( ENVMAP_TYPE_CUBE_UV )
		vec4 texColor = textureCubeUV( envMap, vWorldDirection, backgroundBlurriness );
	#else
		vec4 texColor = vec4( 0.0, 0.0, 0.0, 1.0 );
	#endif
	texColor.rgb *= backgroundIntensity;
	gl_FragColor = texColor;
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
}`, Fu = `varying vec3 vWorldDirection;
#include <common>
void main() {
	vWorldDirection = transformDirection( position, modelMatrix );
	#include <begin_vertex>
	#include <project_vertex>
	gl_Position.z = gl_Position.w;
}`, Bu = `uniform samplerCube tCube;
uniform float tFlip;
uniform float opacity;
varying vec3 vWorldDirection;
void main() {
	vec4 texColor = textureCube( tCube, vec3( tFlip * vWorldDirection.x, vWorldDirection.yz ) );
	gl_FragColor = texColor;
	gl_FragColor.a *= opacity;
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
}`, zu = `#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
varying vec2 vHighPrecisionZW;
void main() {
	#include <uv_vertex>
	#include <skinbase_vertex>
	#ifdef USE_DISPLACEMENTMAP
		#include <beginnormal_vertex>
		#include <morphnormal_vertex>
		#include <skinnormal_vertex>
	#endif
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	vHighPrecisionZW = gl_Position.zw;
}`, Hu = `#if DEPTH_PACKING == 3200
	uniform float opacity;
#endif
#include <common>
#include <packing>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
varying vec2 vHighPrecisionZW;
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( 1.0 );
	#if DEPTH_PACKING == 3200
		diffuseColor.a = opacity;
	#endif
	#include <map_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <logdepthbuf_fragment>
	float fragCoordZ = 0.5 * vHighPrecisionZW[0] / vHighPrecisionZW[1] + 0.5;
	#if DEPTH_PACKING == 3200
		gl_FragColor = vec4( vec3( 1.0 - fragCoordZ ), opacity );
	#elif DEPTH_PACKING == 3201
		gl_FragColor = packDepthToRGBA( fragCoordZ );
	#endif
}`, Gu = `#define DISTANCE
varying vec3 vWorldPosition;
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <skinbase_vertex>
	#ifdef USE_DISPLACEMENTMAP
		#include <beginnormal_vertex>
		#include <morphnormal_vertex>
		#include <skinnormal_vertex>
	#endif
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <worldpos_vertex>
	#include <clipping_planes_vertex>
	vWorldPosition = worldPosition.xyz;
}`, Vu = `#define DISTANCE
uniform vec3 referencePosition;
uniform float nearDistance;
uniform float farDistance;
varying vec3 vWorldPosition;
#include <common>
#include <packing>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <clipping_planes_pars_fragment>
void main () {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( 1.0 );
	#include <map_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	float dist = length( vWorldPosition - referencePosition );
	dist = ( dist - nearDistance ) / ( farDistance - nearDistance );
	dist = saturate( dist );
	gl_FragColor = packDepthToRGBA( dist );
}`, ku = `varying vec3 vWorldDirection;
#include <common>
void main() {
	vWorldDirection = transformDirection( position, modelMatrix );
	#include <begin_vertex>
	#include <project_vertex>
}`, Wu = `uniform sampler2D tEquirect;
varying vec3 vWorldDirection;
#include <common>
void main() {
	vec3 direction = normalize( vWorldDirection );
	vec2 sampleUV = equirectUv( direction );
	gl_FragColor = texture2D( tEquirect, sampleUV );
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
}`, Xu = `uniform float scale;
attribute float lineDistance;
varying float vLineDistance;
#include <common>
#include <uv_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <morphtarget_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	vLineDistance = scale * lineDistance;
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	#include <fog_vertex>
}`, Yu = `uniform vec3 diffuse;
uniform float opacity;
uniform float dashSize;
uniform float totalSize;
varying float vLineDistance;
#include <common>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <fog_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	if ( mod( vLineDistance, totalSize ) > dashSize ) {
		discard;
	}
	vec3 outgoingLight = vec3( 0.0 );
	vec4 diffuseColor = vec4( diffuse, opacity );
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	outgoingLight = diffuseColor.rgb;
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
}`, qu = `#include <common>
#include <uv_pars_vertex>
#include <envmap_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#if defined ( USE_ENVMAP ) || defined ( USE_SKINNING )
		#include <beginnormal_vertex>
		#include <morphnormal_vertex>
		#include <skinbase_vertex>
		#include <skinnormal_vertex>
		#include <defaultnormal_vertex>
	#endif
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	#include <worldpos_vertex>
	#include <envmap_vertex>
	#include <fog_vertex>
}`, ju = `uniform vec3 diffuse;
uniform float opacity;
#ifndef FLAT_SHADED
	varying vec3 vNormal;
#endif
#include <common>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <aomap_pars_fragment>
#include <lightmap_pars_fragment>
#include <envmap_common_pars_fragment>
#include <envmap_pars_fragment>
#include <fog_pars_fragment>
#include <specularmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <specularmap_fragment>
	ReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );
	#ifdef USE_LIGHTMAP
		vec4 lightMapTexel = texture2D( lightMap, vLightMapUv );
		reflectedLight.indirectDiffuse += lightMapTexel.rgb * lightMapIntensity * RECIPROCAL_PI;
	#else
		reflectedLight.indirectDiffuse += vec3( 1.0 );
	#endif
	#include <aomap_fragment>
	reflectedLight.indirectDiffuse *= diffuseColor.rgb;
	vec3 outgoingLight = reflectedLight.indirectDiffuse;
	#include <envmap_fragment>
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, Zu = `#define LAMBERT
varying vec3 vViewPosition;
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <envmap_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	vViewPosition = - mvPosition.xyz;
	#include <worldpos_vertex>
	#include <envmap_vertex>
	#include <shadowmap_vertex>
	#include <fog_vertex>
}`, Ku = `#define LAMBERT
uniform vec3 diffuse;
uniform vec3 emissive;
uniform float opacity;
#include <common>
#include <packing>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <aomap_pars_fragment>
#include <lightmap_pars_fragment>
#include <emissivemap_pars_fragment>
#include <envmap_common_pars_fragment>
#include <envmap_pars_fragment>
#include <fog_pars_fragment>
#include <bsdfs>
#include <lights_pars_begin>
#include <normal_pars_fragment>
#include <lights_lambert_pars_fragment>
#include <shadowmap_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <specularmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	ReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );
	vec3 totalEmissiveRadiance = emissive;
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <specularmap_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	#include <emissivemap_fragment>
	#include <lights_lambert_fragment>
	#include <lights_fragment_begin>
	#include <lights_fragment_maps>
	#include <lights_fragment_end>
	#include <aomap_fragment>
	vec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + totalEmissiveRadiance;
	#include <envmap_fragment>
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, Ju = `#define MATCAP
varying vec3 vViewPosition;
#include <common>
#include <uv_pars_vertex>
#include <color_pars_vertex>
#include <displacementmap_pars_vertex>
#include <fog_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	#include <fog_vertex>
	vViewPosition = - mvPosition.xyz;
}`, $u = `#define MATCAP
uniform vec3 diffuse;
uniform float opacity;
uniform sampler2D matcap;
varying vec3 vViewPosition;
#include <common>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <fog_pars_fragment>
#include <normal_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	vec3 viewDir = normalize( vViewPosition );
	vec3 x = normalize( vec3( viewDir.z, 0.0, - viewDir.x ) );
	vec3 y = cross( viewDir, x );
	vec2 uv = vec2( dot( x, normal ), dot( y, normal ) ) * 0.495 + 0.5;
	#ifdef USE_MATCAP
		vec4 matcapColor = texture2D( matcap, uv );
	#else
		vec4 matcapColor = vec4( vec3( mix( 0.2, 0.8, uv.y ) ), 1.0 );
	#endif
	vec3 outgoingLight = diffuseColor.rgb * matcapColor.rgb;
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, Qu = `#define NORMAL
#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( USE_NORMALMAP_TANGENTSPACE )
	varying vec3 vViewPosition;
#endif
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( USE_NORMALMAP_TANGENTSPACE )
	vViewPosition = - mvPosition.xyz;
#endif
}`, ef = `#define NORMAL
uniform float opacity;
#if defined( FLAT_SHADED ) || defined( USE_BUMPMAP ) || defined( USE_NORMALMAP_TANGENTSPACE )
	varying vec3 vViewPosition;
#endif
#include <packing>
#include <uv_pars_fragment>
#include <normal_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	#include <logdepthbuf_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	gl_FragColor = vec4( packNormalToRGB( normal ), opacity );
	#ifdef OPAQUE
		gl_FragColor.a = 1.0;
	#endif
}`, tf = `#define PHONG
varying vec3 vViewPosition;
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <envmap_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	vViewPosition = - mvPosition.xyz;
	#include <worldpos_vertex>
	#include <envmap_vertex>
	#include <shadowmap_vertex>
	#include <fog_vertex>
}`, nf = `#define PHONG
uniform vec3 diffuse;
uniform vec3 emissive;
uniform vec3 specular;
uniform float shininess;
uniform float opacity;
#include <common>
#include <packing>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <aomap_pars_fragment>
#include <lightmap_pars_fragment>
#include <emissivemap_pars_fragment>
#include <envmap_common_pars_fragment>
#include <envmap_pars_fragment>
#include <fog_pars_fragment>
#include <bsdfs>
#include <lights_pars_begin>
#include <normal_pars_fragment>
#include <lights_phong_pars_fragment>
#include <shadowmap_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <specularmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	ReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );
	vec3 totalEmissiveRadiance = emissive;
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <specularmap_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	#include <emissivemap_fragment>
	#include <lights_phong_fragment>
	#include <lights_fragment_begin>
	#include <lights_fragment_maps>
	#include <lights_fragment_end>
	#include <aomap_fragment>
	vec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + reflectedLight.directSpecular + reflectedLight.indirectSpecular + totalEmissiveRadiance;
	#include <envmap_fragment>
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, rf = `#define STANDARD
varying vec3 vViewPosition;
#ifdef USE_TRANSMISSION
	varying vec3 vWorldPosition;
#endif
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	vViewPosition = - mvPosition.xyz;
	#include <worldpos_vertex>
	#include <shadowmap_vertex>
	#include <fog_vertex>
#ifdef USE_TRANSMISSION
	vWorldPosition = worldPosition.xyz;
#endif
}`, sf = `#define STANDARD
#ifdef PHYSICAL
	#define IOR
	#define USE_SPECULAR
#endif
uniform vec3 diffuse;
uniform vec3 emissive;
uniform float roughness;
uniform float metalness;
uniform float opacity;
#ifdef IOR
	uniform float ior;
#endif
#ifdef USE_SPECULAR
	uniform float specularIntensity;
	uniform vec3 specularColor;
	#ifdef USE_SPECULAR_COLORMAP
		uniform sampler2D specularColorMap;
	#endif
	#ifdef USE_SPECULAR_INTENSITYMAP
		uniform sampler2D specularIntensityMap;
	#endif
#endif
#ifdef USE_CLEARCOAT
	uniform float clearcoat;
	uniform float clearcoatRoughness;
#endif
#ifdef USE_IRIDESCENCE
	uniform float iridescence;
	uniform float iridescenceIOR;
	uniform float iridescenceThicknessMinimum;
	uniform float iridescenceThicknessMaximum;
#endif
#ifdef USE_SHEEN
	uniform vec3 sheenColor;
	uniform float sheenRoughness;
	#ifdef USE_SHEEN_COLORMAP
		uniform sampler2D sheenColorMap;
	#endif
	#ifdef USE_SHEEN_ROUGHNESSMAP
		uniform sampler2D sheenRoughnessMap;
	#endif
#endif
#ifdef USE_ANISOTROPY
	uniform vec2 anisotropyVector;
	#ifdef USE_ANISOTROPYMAP
		uniform sampler2D anisotropyMap;
	#endif
#endif
varying vec3 vViewPosition;
#include <common>
#include <packing>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <aomap_pars_fragment>
#include <lightmap_pars_fragment>
#include <emissivemap_pars_fragment>
#include <iridescence_fragment>
#include <cube_uv_reflection_fragment>
#include <envmap_common_pars_fragment>
#include <envmap_physical_pars_fragment>
#include <fog_pars_fragment>
#include <lights_pars_begin>
#include <normal_pars_fragment>
#include <lights_physical_pars_fragment>
#include <transmission_pars_fragment>
#include <shadowmap_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <clearcoat_pars_fragment>
#include <iridescence_pars_fragment>
#include <roughnessmap_pars_fragment>
#include <metalnessmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	ReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );
	vec3 totalEmissiveRadiance = emissive;
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <roughnessmap_fragment>
	#include <metalnessmap_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	#include <clearcoat_normal_fragment_begin>
	#include <clearcoat_normal_fragment_maps>
	#include <emissivemap_fragment>
	#include <lights_physical_fragment>
	#include <lights_fragment_begin>
	#include <lights_fragment_maps>
	#include <lights_fragment_end>
	#include <aomap_fragment>
	vec3 totalDiffuse = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse;
	vec3 totalSpecular = reflectedLight.directSpecular + reflectedLight.indirectSpecular;
	#include <transmission_fragment>
	vec3 outgoingLight = totalDiffuse + totalSpecular + totalEmissiveRadiance;
	#ifdef USE_SHEEN
		float sheenEnergyComp = 1.0 - 0.157 * max3( material.sheenColor );
		outgoingLight = outgoingLight * sheenEnergyComp + sheenSpecular;
	#endif
	#ifdef USE_CLEARCOAT
		float dotNVcc = saturate( dot( geometry.clearcoatNormal, geometry.viewDir ) );
		vec3 Fcc = F_Schlick( material.clearcoatF0, material.clearcoatF90, dotNVcc );
		outgoingLight = outgoingLight * ( 1.0 - material.clearcoat * Fcc ) + clearcoatSpecular * material.clearcoat;
	#endif
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, af = `#define TOON
varying vec3 vViewPosition;
#include <common>
#include <uv_pars_vertex>
#include <displacementmap_pars_vertex>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <normal_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <normal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <displacementmap_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	vViewPosition = - mvPosition.xyz;
	#include <worldpos_vertex>
	#include <shadowmap_vertex>
	#include <fog_vertex>
}`, of = `#define TOON
uniform vec3 diffuse;
uniform vec3 emissive;
uniform float opacity;
#include <common>
#include <packing>
#include <dithering_pars_fragment>
#include <color_pars_fragment>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <aomap_pars_fragment>
#include <lightmap_pars_fragment>
#include <emissivemap_pars_fragment>
#include <gradientmap_pars_fragment>
#include <fog_pars_fragment>
#include <bsdfs>
#include <lights_pars_begin>
#include <normal_pars_fragment>
#include <lights_toon_pars_fragment>
#include <shadowmap_pars_fragment>
#include <bumpmap_pars_fragment>
#include <normalmap_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec4 diffuseColor = vec4( diffuse, opacity );
	ReflectedLight reflectedLight = ReflectedLight( vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ), vec3( 0.0 ) );
	vec3 totalEmissiveRadiance = emissive;
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <color_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	#include <normal_fragment_begin>
	#include <normal_fragment_maps>
	#include <emissivemap_fragment>
	#include <lights_toon_fragment>
	#include <lights_fragment_begin>
	#include <lights_fragment_maps>
	#include <lights_fragment_end>
	#include <aomap_fragment>
	vec3 outgoingLight = reflectedLight.directDiffuse + reflectedLight.indirectDiffuse + totalEmissiveRadiance;
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
	#include <dithering_fragment>
}`, lf = `uniform float size;
uniform float scale;
#include <common>
#include <color_pars_vertex>
#include <fog_pars_vertex>
#include <morphtarget_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
#ifdef USE_POINTS_UV
	varying vec2 vUv;
	uniform mat3 uvTransform;
#endif
void main() {
	#ifdef USE_POINTS_UV
		vUv = ( uvTransform * vec3( uv, 1 ) ).xy;
	#endif
	#include <color_vertex>
	#include <morphcolor_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <project_vertex>
	gl_PointSize = size;
	#ifdef USE_SIZEATTENUATION
		bool isPerspective = isPerspectiveMatrix( projectionMatrix );
		if ( isPerspective ) gl_PointSize *= ( scale / - mvPosition.z );
	#endif
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	#include <worldpos_vertex>
	#include <fog_vertex>
}`, cf = `uniform vec3 diffuse;
uniform float opacity;
#include <common>
#include <color_pars_fragment>
#include <map_particle_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <fog_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec3 outgoingLight = vec3( 0.0 );
	vec4 diffuseColor = vec4( diffuse, opacity );
	#include <logdepthbuf_fragment>
	#include <map_particle_fragment>
	#include <color_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	outgoingLight = diffuseColor.rgb;
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
	#include <premultiplied_alpha_fragment>
}`, hf = `#include <common>
#include <fog_pars_vertex>
#include <morphtarget_pars_vertex>
#include <skinning_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <shadowmap_pars_vertex>
void main() {
	#include <beginnormal_vertex>
	#include <morphnormal_vertex>
	#include <skinbase_vertex>
	#include <skinnormal_vertex>
	#include <defaultnormal_vertex>
	#include <begin_vertex>
	#include <morphtarget_vertex>
	#include <skinning_vertex>
	#include <project_vertex>
	#include <logdepthbuf_vertex>
	#include <worldpos_vertex>
	#include <shadowmap_vertex>
	#include <fog_vertex>
}`, uf = `uniform vec3 color;
uniform float opacity;
#include <common>
#include <packing>
#include <fog_pars_fragment>
#include <bsdfs>
#include <lights_pars_begin>
#include <logdepthbuf_pars_fragment>
#include <shadowmap_pars_fragment>
#include <shadowmask_pars_fragment>
void main() {
	#include <logdepthbuf_fragment>
	gl_FragColor = vec4( color, opacity * ( 1.0 - getShadowMask() ) );
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
}`, ff = `uniform float rotation;
uniform vec2 center;
#include <common>
#include <uv_pars_vertex>
#include <fog_pars_vertex>
#include <logdepthbuf_pars_vertex>
#include <clipping_planes_pars_vertex>
void main() {
	#include <uv_vertex>
	vec4 mvPosition = modelViewMatrix * vec4( 0.0, 0.0, 0.0, 1.0 );
	vec2 scale;
	scale.x = length( vec3( modelMatrix[ 0 ].x, modelMatrix[ 0 ].y, modelMatrix[ 0 ].z ) );
	scale.y = length( vec3( modelMatrix[ 1 ].x, modelMatrix[ 1 ].y, modelMatrix[ 1 ].z ) );
	#ifndef USE_SIZEATTENUATION
		bool isPerspective = isPerspectiveMatrix( projectionMatrix );
		if ( isPerspective ) scale *= - mvPosition.z;
	#endif
	vec2 alignedPosition = ( position.xy - ( center - vec2( 0.5 ) ) ) * scale;
	vec2 rotatedPosition;
	rotatedPosition.x = cos( rotation ) * alignedPosition.x - sin( rotation ) * alignedPosition.y;
	rotatedPosition.y = sin( rotation ) * alignedPosition.x + cos( rotation ) * alignedPosition.y;
	mvPosition.xy += rotatedPosition;
	gl_Position = projectionMatrix * mvPosition;
	#include <logdepthbuf_vertex>
	#include <clipping_planes_vertex>
	#include <fog_vertex>
}`, df = `uniform vec3 diffuse;
uniform float opacity;
#include <common>
#include <uv_pars_fragment>
#include <map_pars_fragment>
#include <alphamap_pars_fragment>
#include <alphatest_pars_fragment>
#include <alphahash_pars_fragment>
#include <fog_pars_fragment>
#include <logdepthbuf_pars_fragment>
#include <clipping_planes_pars_fragment>
void main() {
	#include <clipping_planes_fragment>
	vec3 outgoingLight = vec3( 0.0 );
	vec4 diffuseColor = vec4( diffuse, opacity );
	#include <logdepthbuf_fragment>
	#include <map_fragment>
	#include <alphamap_fragment>
	#include <alphatest_fragment>
	#include <alphahash_fragment>
	outgoingLight = diffuseColor.rgb;
	#include <opaque_fragment>
	#include <tonemapping_fragment>
	#include <colorspace_fragment>
	#include <fog_fragment>
}`, Ne = {
  alphahash_fragment: Oc,
  alphahash_pars_fragment: Fc,
  alphamap_fragment: Bc,
  alphamap_pars_fragment: zc,
  alphatest_fragment: Hc,
  alphatest_pars_fragment: Gc,
  aomap_fragment: Vc,
  aomap_pars_fragment: kc,
  begin_vertex: Wc,
  beginnormal_vertex: Xc,
  bsdfs: Yc,
  iridescence_fragment: qc,
  bumpmap_pars_fragment: jc,
  clipping_planes_fragment: Zc,
  clipping_planes_pars_fragment: Kc,
  clipping_planes_pars_vertex: Jc,
  clipping_planes_vertex: $c,
  color_fragment: Qc,
  color_pars_fragment: eh,
  color_pars_vertex: th,
  color_vertex: nh,
  common: ih,
  cube_uv_reflection_fragment: rh,
  defaultnormal_vertex: sh,
  displacementmap_pars_vertex: ah,
  displacementmap_vertex: oh,
  emissivemap_fragment: lh,
  emissivemap_pars_fragment: ch,
  colorspace_fragment: hh,
  colorspace_pars_fragment: uh,
  envmap_fragment: fh,
  envmap_common_pars_fragment: dh,
  envmap_pars_fragment: ph,
  envmap_pars_vertex: mh,
  envmap_physical_pars_fragment: wh,
  envmap_vertex: gh,
  fog_vertex: _h,
  fog_pars_vertex: vh,
  fog_fragment: xh,
  fog_pars_fragment: Mh,
  gradientmap_pars_fragment: Sh,
  lightmap_fragment: Eh,
  lightmap_pars_fragment: yh,
  lights_lambert_fragment: Th,
  lights_lambert_pars_fragment: Ah,
  lights_pars_begin: bh,
  lights_toon_fragment: Rh,
  lights_toon_pars_fragment: Ch,
  lights_phong_fragment: Ph,
  lights_phong_pars_fragment: Lh,
  lights_physical_fragment: Uh,
  lights_physical_pars_fragment: Dh,
  lights_fragment_begin: Ih,
  lights_fragment_maps: Nh,
  lights_fragment_end: Oh,
  logdepthbuf_fragment: Fh,
  logdepthbuf_pars_fragment: Bh,
  logdepthbuf_pars_vertex: zh,
  logdepthbuf_vertex: Hh,
  map_fragment: Gh,
  map_pars_fragment: Vh,
  map_particle_fragment: kh,
  map_particle_pars_fragment: Wh,
  metalnessmap_fragment: Xh,
  metalnessmap_pars_fragment: Yh,
  morphcolor_vertex: qh,
  morphnormal_vertex: jh,
  morphtarget_pars_vertex: Zh,
  morphtarget_vertex: Kh,
  normal_fragment_begin: Jh,
  normal_fragment_maps: $h,
  normal_pars_fragment: Qh,
  normal_pars_vertex: eu,
  normal_vertex: tu,
  normalmap_pars_fragment: nu,
  clearcoat_normal_fragment_begin: iu,
  clearcoat_normal_fragment_maps: ru,
  clearcoat_pars_fragment: su,
  iridescence_pars_fragment: au,
  opaque_fragment: ou,
  packing: lu,
  premultiplied_alpha_fragment: cu,
  project_vertex: hu,
  dithering_fragment: uu,
  dithering_pars_fragment: fu,
  roughnessmap_fragment: du,
  roughnessmap_pars_fragment: pu,
  shadowmap_pars_fragment: mu,
  shadowmap_pars_vertex: gu,
  shadowmap_vertex: _u,
  shadowmask_pars_fragment: vu,
  skinbase_vertex: xu,
  skinning_pars_vertex: Mu,
  skinning_vertex: Su,
  skinnormal_vertex: Eu,
  specularmap_fragment: yu,
  specularmap_pars_fragment: Tu,
  tonemapping_fragment: Au,
  tonemapping_pars_fragment: bu,
  transmission_fragment: wu,
  transmission_pars_fragment: Ru,
  uv_pars_fragment: Cu,
  uv_pars_vertex: Pu,
  uv_vertex: Lu,
  worldpos_vertex: Uu,
  background_vert: Du,
  background_frag: Iu,
  backgroundCube_vert: Nu,
  backgroundCube_frag: Ou,
  cube_vert: Fu,
  cube_frag: Bu,
  depth_vert: zu,
  depth_frag: Hu,
  distanceRGBA_vert: Gu,
  distanceRGBA_frag: Vu,
  equirect_vert: ku,
  equirect_frag: Wu,
  linedashed_vert: Xu,
  linedashed_frag: Yu,
  meshbasic_vert: qu,
  meshbasic_frag: ju,
  meshlambert_vert: Zu,
  meshlambert_frag: Ku,
  meshmatcap_vert: Ju,
  meshmatcap_frag: $u,
  meshnormal_vert: Qu,
  meshnormal_frag: ef,
  meshphong_vert: tf,
  meshphong_frag: nf,
  meshphysical_vert: rf,
  meshphysical_frag: sf,
  meshtoon_vert: af,
  meshtoon_frag: of,
  points_vert: lf,
  points_frag: cf,
  shadow_vert: hf,
  shadow_frag: uf,
  sprite_vert: ff,
  sprite_frag: df
}, ue = {
  common: {
    diffuse: { value: /* @__PURE__ */ new We(16777215) },
    opacity: { value: 1 },
    map: { value: null },
    mapTransform: { value: /* @__PURE__ */ new Be() },
    alphaMap: { value: null },
    alphaMapTransform: { value: /* @__PURE__ */ new Be() },
    alphaTest: { value: 0 }
  },
  specularmap: {
    specularMap: { value: null },
    specularMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  envmap: {
    envMap: { value: null },
    flipEnvMap: { value: -1 },
    reflectivity: { value: 1 },
    // basic, lambert, phong
    ior: { value: 1.5 },
    // physical
    refractionRatio: { value: 0.98 }
    // basic, lambert, phong
  },
  aomap: {
    aoMap: { value: null },
    aoMapIntensity: { value: 1 },
    aoMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  lightmap: {
    lightMap: { value: null },
    lightMapIntensity: { value: 1 },
    lightMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  bumpmap: {
    bumpMap: { value: null },
    bumpMapTransform: { value: /* @__PURE__ */ new Be() },
    bumpScale: { value: 1 }
  },
  normalmap: {
    normalMap: { value: null },
    normalMapTransform: { value: /* @__PURE__ */ new Be() },
    normalScale: { value: /* @__PURE__ */ new oe(1, 1) }
  },
  displacementmap: {
    displacementMap: { value: null },
    displacementMapTransform: { value: /* @__PURE__ */ new Be() },
    displacementScale: { value: 1 },
    displacementBias: { value: 0 }
  },
  emissivemap: {
    emissiveMap: { value: null },
    emissiveMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  metalnessmap: {
    metalnessMap: { value: null },
    metalnessMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  roughnessmap: {
    roughnessMap: { value: null },
    roughnessMapTransform: { value: /* @__PURE__ */ new Be() }
  },
  gradientmap: {
    gradientMap: { value: null }
  },
  fog: {
    fogDensity: { value: 25e-5 },
    fogNear: { value: 1 },
    fogFar: { value: 2e3 },
    fogColor: { value: /* @__PURE__ */ new We(16777215) }
  },
  lights: {
    ambientLightColor: { value: [] },
    lightProbe: { value: [] },
    directionalLights: { value: [], properties: {
      direction: {},
      color: {}
    } },
    directionalLightShadows: { value: [], properties: {
      shadowBias: {},
      shadowNormalBias: {},
      shadowRadius: {},
      shadowMapSize: {}
    } },
    directionalShadowMap: { value: [] },
    directionalShadowMatrix: { value: [] },
    spotLights: { value: [], properties: {
      color: {},
      position: {},
      direction: {},
      distance: {},
      coneCos: {},
      penumbraCos: {},
      decay: {}
    } },
    spotLightShadows: { value: [], properties: {
      shadowBias: {},
      shadowNormalBias: {},
      shadowRadius: {},
      shadowMapSize: {}
    } },
    spotLightMap: { value: [] },
    spotShadowMap: { value: [] },
    spotLightMatrix: { value: [] },
    pointLights: { value: [], properties: {
      color: {},
      position: {},
      decay: {},
      distance: {}
    } },
    pointLightShadows: { value: [], properties: {
      shadowBias: {},
      shadowNormalBias: {},
      shadowRadius: {},
      shadowMapSize: {},
      shadowCameraNear: {},
      shadowCameraFar: {}
    } },
    pointShadowMap: { value: [] },
    pointShadowMatrix: { value: [] },
    hemisphereLights: { value: [], properties: {
      direction: {},
      skyColor: {},
      groundColor: {}
    } },
    // TODO (abelnation): RectAreaLight BRDF data needs to be moved from example to main src
    rectAreaLights: { value: [], properties: {
      color: {},
      position: {},
      width: {},
      height: {}
    } },
    ltc_1: { value: null },
    ltc_2: { value: null }
  },
  points: {
    diffuse: { value: /* @__PURE__ */ new We(16777215) },
    opacity: { value: 1 },
    size: { value: 1 },
    scale: { value: 1 },
    map: { value: null },
    alphaMap: { value: null },
    alphaMapTransform: { value: /* @__PURE__ */ new Be() },
    alphaTest: { value: 0 },
    uvTransform: { value: /* @__PURE__ */ new Be() }
  },
  sprite: {
    diffuse: { value: /* @__PURE__ */ new We(16777215) },
    opacity: { value: 1 },
    center: { value: /* @__PURE__ */ new oe(0.5, 0.5) },
    rotation: { value: 0 },
    map: { value: null },
    mapTransform: { value: /* @__PURE__ */ new Be() },
    alphaMap: { value: null },
    alphaMapTransform: { value: /* @__PURE__ */ new Be() },
    alphaTest: { value: 0 }
  }
}, It = {
  basic: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.specularmap,
      ue.envmap,
      ue.aomap,
      ue.lightmap,
      ue.fog
    ]),
    vertexShader: Ne.meshbasic_vert,
    fragmentShader: Ne.meshbasic_frag
  },
  lambert: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.specularmap,
      ue.envmap,
      ue.aomap,
      ue.lightmap,
      ue.emissivemap,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      ue.fog,
      ue.lights,
      {
        emissive: { value: /* @__PURE__ */ new We(0) }
      }
    ]),
    vertexShader: Ne.meshlambert_vert,
    fragmentShader: Ne.meshlambert_frag
  },
  phong: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.specularmap,
      ue.envmap,
      ue.aomap,
      ue.lightmap,
      ue.emissivemap,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      ue.fog,
      ue.lights,
      {
        emissive: { value: /* @__PURE__ */ new We(0) },
        specular: { value: /* @__PURE__ */ new We(1118481) },
        shininess: { value: 30 }
      }
    ]),
    vertexShader: Ne.meshphong_vert,
    fragmentShader: Ne.meshphong_frag
  },
  standard: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.envmap,
      ue.aomap,
      ue.lightmap,
      ue.emissivemap,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      ue.roughnessmap,
      ue.metalnessmap,
      ue.fog,
      ue.lights,
      {
        emissive: { value: /* @__PURE__ */ new We(0) },
        roughness: { value: 1 },
        metalness: { value: 0 },
        envMapIntensity: { value: 1 }
        // temporary
      }
    ]),
    vertexShader: Ne.meshphysical_vert,
    fragmentShader: Ne.meshphysical_frag
  },
  toon: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.aomap,
      ue.lightmap,
      ue.emissivemap,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      ue.gradientmap,
      ue.fog,
      ue.lights,
      {
        emissive: { value: /* @__PURE__ */ new We(0) }
      }
    ]),
    vertexShader: Ne.meshtoon_vert,
    fragmentShader: Ne.meshtoon_frag
  },
  matcap: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      ue.fog,
      {
        matcap: { value: null }
      }
    ]),
    vertexShader: Ne.meshmatcap_vert,
    fragmentShader: Ne.meshmatcap_frag
  },
  points: {
    uniforms: /* @__PURE__ */ pt([
      ue.points,
      ue.fog
    ]),
    vertexShader: Ne.points_vert,
    fragmentShader: Ne.points_frag
  },
  dashed: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.fog,
      {
        scale: { value: 1 },
        dashSize: { value: 1 },
        totalSize: { value: 2 }
      }
    ]),
    vertexShader: Ne.linedashed_vert,
    fragmentShader: Ne.linedashed_frag
  },
  depth: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.displacementmap
    ]),
    vertexShader: Ne.depth_vert,
    fragmentShader: Ne.depth_frag
  },
  normal: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.bumpmap,
      ue.normalmap,
      ue.displacementmap,
      {
        opacity: { value: 1 }
      }
    ]),
    vertexShader: Ne.meshnormal_vert,
    fragmentShader: Ne.meshnormal_frag
  },
  sprite: {
    uniforms: /* @__PURE__ */ pt([
      ue.sprite,
      ue.fog
    ]),
    vertexShader: Ne.sprite_vert,
    fragmentShader: Ne.sprite_frag
  },
  background: {
    uniforms: {
      uvTransform: { value: /* @__PURE__ */ new Be() },
      t2D: { value: null },
      backgroundIntensity: { value: 1 }
    },
    vertexShader: Ne.background_vert,
    fragmentShader: Ne.background_frag
  },
  backgroundCube: {
    uniforms: {
      envMap: { value: null },
      flipEnvMap: { value: -1 },
      backgroundBlurriness: { value: 0 },
      backgroundIntensity: { value: 1 }
    },
    vertexShader: Ne.backgroundCube_vert,
    fragmentShader: Ne.backgroundCube_frag
  },
  cube: {
    uniforms: {
      tCube: { value: null },
      tFlip: { value: -1 },
      opacity: { value: 1 }
    },
    vertexShader: Ne.cube_vert,
    fragmentShader: Ne.cube_frag
  },
  equirect: {
    uniforms: {
      tEquirect: { value: null }
    },
    vertexShader: Ne.equirect_vert,
    fragmentShader: Ne.equirect_frag
  },
  distanceRGBA: {
    uniforms: /* @__PURE__ */ pt([
      ue.common,
      ue.displacementmap,
      {
        referencePosition: { value: /* @__PURE__ */ new U() },
        nearDistance: { value: 1 },
        farDistance: { value: 1e3 }
      }
    ]),
    vertexShader: Ne.distanceRGBA_vert,
    fragmentShader: Ne.distanceRGBA_frag
  },
  shadow: {
    uniforms: /* @__PURE__ */ pt([
      ue.lights,
      ue.fog,
      {
        color: { value: /* @__PURE__ */ new We(0) },
        opacity: { value: 1 }
      }
    ]),
    vertexShader: Ne.shadow_vert,
    fragmentShader: Ne.shadow_frag
  }
};
It.physical = {
  uniforms: /* @__PURE__ */ pt([
    It.standard.uniforms,
    {
      clearcoat: { value: 0 },
      clearcoatMap: { value: null },
      clearcoatMapTransform: { value: /* @__PURE__ */ new Be() },
      clearcoatNormalMap: { value: null },
      clearcoatNormalMapTransform: { value: /* @__PURE__ */ new Be() },
      clearcoatNormalScale: { value: /* @__PURE__ */ new oe(1, 1) },
      clearcoatRoughness: { value: 0 },
      clearcoatRoughnessMap: { value: null },
      clearcoatRoughnessMapTransform: { value: /* @__PURE__ */ new Be() },
      iridescence: { value: 0 },
      iridescenceMap: { value: null },
      iridescenceMapTransform: { value: /* @__PURE__ */ new Be() },
      iridescenceIOR: { value: 1.3 },
      iridescenceThicknessMinimum: { value: 100 },
      iridescenceThicknessMaximum: { value: 400 },
      iridescenceThicknessMap: { value: null },
      iridescenceThicknessMapTransform: { value: /* @__PURE__ */ new Be() },
      sheen: { value: 0 },
      sheenColor: { value: /* @__PURE__ */ new We(0) },
      sheenColorMap: { value: null },
      sheenColorMapTransform: { value: /* @__PURE__ */ new Be() },
      sheenRoughness: { value: 1 },
      sheenRoughnessMap: { value: null },
      sheenRoughnessMapTransform: { value: /* @__PURE__ */ new Be() },
      transmission: { value: 0 },
      transmissionMap: { value: null },
      transmissionMapTransform: { value: /* @__PURE__ */ new Be() },
      transmissionSamplerSize: { value: /* @__PURE__ */ new oe() },
      transmissionSamplerMap: { value: null },
      thickness: { value: 0 },
      thicknessMap: { value: null },
      thicknessMapTransform: { value: /* @__PURE__ */ new Be() },
      attenuationDistance: { value: 0 },
      attenuationColor: { value: /* @__PURE__ */ new We(0) },
      specularColor: { value: /* @__PURE__ */ new We(1, 1, 1) },
      specularColorMap: { value: null },
      specularColorMapTransform: { value: /* @__PURE__ */ new Be() },
      specularIntensity: { value: 1 },
      specularIntensityMap: { value: null },
      specularIntensityMapTransform: { value: /* @__PURE__ */ new Be() },
      anisotropyVector: { value: /* @__PURE__ */ new oe() },
      anisotropyMap: { value: null },
      anisotropyMapTransform: { value: /* @__PURE__ */ new Be() }
    }
  ]),
  vertexShader: Ne.meshphysical_vert,
  fragmentShader: Ne.meshphysical_frag
};
const qi = { r: 0, b: 0, g: 0 };
function pf(i, e, t, n, r, s, o) {
  const a = new We(0);
  let l = s === !0 ? 0 : 1, c, h, f = null, u = 0, m = null;
  function g(p, d) {
    let A = !1, _ = d.isScene === !0 ? d.background : null;
    switch (_ && _.isTexture && (_ = (d.backgroundBlurriness > 0 ? t : e).get(_)), _ === null ? x(a, l) : _ && _.isColor && (x(_, 1), A = !0), i.xr.getEnvironmentBlendMode()) {
      case "opaque":
        A = !0;
        break;
      case "additive":
        n.buffers.color.setClear(0, 0, 0, 1, o), A = !0;
        break;
      case "alpha-blend":
        n.buffers.color.setClear(0, 0, 0, 0, o), A = !0;
        break;
    }
    (i.autoClear || A) && i.clear(i.autoClearColor, i.autoClearDepth, i.autoClearStencil), _ && (_.isCubeTexture || _.mapping === rr) ? (h === void 0 && (h = new Nt(
      new bi(1, 1, 1),
      new bn({
        name: "BackgroundCubeMaterial",
        uniforms: ii(It.backgroundCube.uniforms),
        vertexShader: It.backgroundCube.vertexShader,
        fragmentShader: It.backgroundCube.fragmentShader,
        side: gt,
        depthTest: !1,
        depthWrite: !1,
        fog: !1
      })
    ), h.geometry.deleteAttribute("normal"), h.geometry.deleteAttribute("uv"), h.onBeforeRender = function(L, w, V) {
      this.matrixWorld.copyPosition(V.matrixWorld);
    }, Object.defineProperty(h.material, "envMap", {
      get: function() {
        return this.uniforms.envMap.value;
      }
    }), r.update(h)), h.material.uniforms.envMap.value = _, h.material.uniforms.flipEnvMap.value = _.isCubeTexture && _.isRenderTargetTexture === !1 ? -1 : 1, h.material.uniforms.backgroundBlurriness.value = d.backgroundBlurriness, h.material.uniforms.backgroundIntensity.value = d.backgroundIntensity, h.material.toneMapped = _.colorSpace !== Oe, (f !== _ || u !== _.version || m !== i.toneMapping) && (h.material.needsUpdate = !0, f = _, u = _.version, m = i.toneMapping), h.layers.enableAll(), p.unshift(h, h.geometry, h.material, 0, 0, null)) : _ && _.isTexture && (c === void 0 && (c = new Nt(
      new ar(2, 2),
      new bn({
        name: "BackgroundMaterial",
        uniforms: ii(It.background.uniforms),
        vertexShader: It.background.vertexShader,
        fragmentShader: It.background.fragmentShader,
        side: ln,
        depthTest: !1,
        depthWrite: !1,
        fog: !1
      })
    ), c.geometry.deleteAttribute("normal"), Object.defineProperty(c.material, "map", {
      get: function() {
        return this.uniforms.t2D.value;
      }
    }), r.update(c)), c.material.uniforms.t2D.value = _, c.material.uniforms.backgroundIntensity.value = d.backgroundIntensity, c.material.toneMapped = _.colorSpace !== Oe, _.matrixAutoUpdate === !0 && _.updateMatrix(), c.material.uniforms.uvTransform.value.copy(_.matrix), (f !== _ || u !== _.version || m !== i.toneMapping) && (c.material.needsUpdate = !0, f = _, u = _.version, m = i.toneMapping), c.layers.enableAll(), p.unshift(c, c.geometry, c.material, 0, 0, null));
  }
  function x(p, d) {
    p.getRGB(qi, So(i)), n.buffers.color.setClear(qi.r, qi.g, qi.b, d, o);
  }
  return {
    getClearColor: function() {
      return a;
    },
    setClearColor: function(p, d = 1) {
      a.set(p), l = d, x(a, l);
    },
    getClearAlpha: function() {
      return l;
    },
    setClearAlpha: function(p) {
      l = p, x(a, l);
    },
    render: g
  };
}
function mf(i, e, t, n) {
  const r = i.getParameter(i.MAX_VERTEX_ATTRIBS), s = n.isWebGL2 ? null : e.get("OES_vertex_array_object"), o = n.isWebGL2 || s !== null, a = {}, l = p(null);
  let c = l, h = !1;
  function f(H, G, Q, X, j) {
    let J = !1;
    if (o) {
      const ee = x(X, Q, G);
      c !== ee && (c = ee, m(c.object)), J = d(H, X, Q, j), J && A(H, X, Q, j);
    } else {
      const ee = G.wireframe === !0;
      (c.geometry !== X.id || c.program !== Q.id || c.wireframe !== ee) && (c.geometry = X.id, c.program = Q.id, c.wireframe = ee, J = !0);
    }
    j !== null && t.update(j, i.ELEMENT_ARRAY_BUFFER), (J || h) && (h = !1, V(H, G, Q, X), j !== null && i.bindBuffer(i.ELEMENT_ARRAY_BUFFER, t.get(j).buffer));
  }
  function u() {
    return n.isWebGL2 ? i.createVertexArray() : s.createVertexArrayOES();
  }
  function m(H) {
    return n.isWebGL2 ? i.bindVertexArray(H) : s.bindVertexArrayOES(H);
  }
  function g(H) {
    return n.isWebGL2 ? i.deleteVertexArray(H) : s.deleteVertexArrayOES(H);
  }
  function x(H, G, Q) {
    const X = Q.wireframe === !0;
    let j = a[H.id];
    j === void 0 && (j = {}, a[H.id] = j);
    let J = j[G.id];
    J === void 0 && (J = {}, j[G.id] = J);
    let ee = J[X];
    return ee === void 0 && (ee = p(u()), J[X] = ee), ee;
  }
  function p(H) {
    const G = [], Q = [], X = [];
    for (let j = 0; j < r; j++)
      G[j] = 0, Q[j] = 0, X[j] = 0;
    return {
      // for backward compatibility on non-VAO support browser
      geometry: null,
      program: null,
      wireframe: !1,
      newAttributes: G,
      enabledAttributes: Q,
      attributeDivisors: X,
      object: H,
      attributes: {},
      index: null
    };
  }
  function d(H, G, Q, X) {
    const j = c.attributes, J = G.attributes;
    let ee = 0;
    const I = Q.getAttributes();
    for (const q in I)
      if (I[q].location >= 0) {
        const me = j[q];
        let ve = J[q];
        if (ve === void 0 && (q === "instanceMatrix" && H.instanceMatrix && (ve = H.instanceMatrix), q === "instanceColor" && H.instanceColor && (ve = H.instanceColor)), me === void 0 || me.attribute !== ve || ve && me.data !== ve.data)
          return !0;
        ee++;
      }
    return c.attributesNum !== ee || c.index !== X;
  }
  function A(H, G, Q, X) {
    const j = {}, J = G.attributes;
    let ee = 0;
    const I = Q.getAttributes();
    for (const q in I)
      if (I[q].location >= 0) {
        let me = J[q];
        me === void 0 && (q === "instanceMatrix" && H.instanceMatrix && (me = H.instanceMatrix), q === "instanceColor" && H.instanceColor && (me = H.instanceColor));
        const ve = {};
        ve.attribute = me, me && me.data && (ve.data = me.data), j[q] = ve, ee++;
      }
    c.attributes = j, c.attributesNum = ee, c.index = X;
  }
  function _() {
    const H = c.newAttributes;
    for (let G = 0, Q = H.length; G < Q; G++)
      H[G] = 0;
  }
  function T(H) {
    C(H, 0);
  }
  function C(H, G) {
    const Q = c.newAttributes, X = c.enabledAttributes, j = c.attributeDivisors;
    Q[H] = 1, X[H] === 0 && (i.enableVertexAttribArray(H), X[H] = 1), j[H] !== G && ((n.isWebGL2 ? i : e.get("ANGLE_instanced_arrays"))[n.isWebGL2 ? "vertexAttribDivisor" : "vertexAttribDivisorANGLE"](H, G), j[H] = G);
  }
  function L() {
    const H = c.newAttributes, G = c.enabledAttributes;
    for (let Q = 0, X = G.length; Q < X; Q++)
      G[Q] !== H[Q] && (i.disableVertexAttribArray(Q), G[Q] = 0);
  }
  function w(H, G, Q, X, j, J, ee) {
    ee === !0 ? i.vertexAttribIPointer(H, G, Q, j, J) : i.vertexAttribPointer(H, G, Q, X, j, J);
  }
  function V(H, G, Q, X) {
    if (n.isWebGL2 === !1 && (H.isInstancedMesh || X.isInstancedBufferGeometry) && e.get("ANGLE_instanced_arrays") === null)
      return;
    _();
    const j = X.attributes, J = Q.getAttributes(), ee = G.defaultAttributeValues;
    for (const I in J) {
      const q = J[I];
      if (q.location >= 0) {
        let pe = j[I];
        if (pe === void 0 && (I === "instanceMatrix" && H.instanceMatrix && (pe = H.instanceMatrix), I === "instanceColor" && H.instanceColor && (pe = H.instanceColor)), pe !== void 0) {
          const me = pe.normalized, ve = pe.itemSize, Ae = t.get(pe);
          if (Ae === void 0)
            continue;
          const be = Ae.buffer, Te = Ae.type, ke = Ae.bytesPerElement, Ye = n.isWebGL2 === !0 && (Te === i.INT || Te === i.UNSIGNED_INT || pe.gpuType === io);
          if (pe.isInterleavedBufferAttribute) {
            const we = pe.data, b = we.stride, le = pe.offset;
            if (we.isInstancedInterleavedBuffer) {
              for (let Z = 0; Z < q.locationSize; Z++)
                C(q.location + Z, we.meshPerAttribute);
              H.isInstancedMesh !== !0 && X._maxInstanceCount === void 0 && (X._maxInstanceCount = we.meshPerAttribute * we.count);
            } else
              for (let Z = 0; Z < q.locationSize; Z++)
                T(q.location + Z);
            i.bindBuffer(i.ARRAY_BUFFER, be);
            for (let Z = 0; Z < q.locationSize; Z++)
              w(
                q.location + Z,
                ve / q.locationSize,
                Te,
                me,
                b * ke,
                (le + ve / q.locationSize * Z) * ke,
                Ye
              );
          } else {
            if (pe.isInstancedBufferAttribute) {
              for (let we = 0; we < q.locationSize; we++)
                C(q.location + we, pe.meshPerAttribute);
              H.isInstancedMesh !== !0 && X._maxInstanceCount === void 0 && (X._maxInstanceCount = pe.meshPerAttribute * pe.count);
            } else
              for (let we = 0; we < q.locationSize; we++)
                T(q.location + we);
            i.bindBuffer(i.ARRAY_BUFFER, be);
            for (let we = 0; we < q.locationSize; we++)
              w(
                q.location + we,
                ve / q.locationSize,
                Te,
                me,
                ve * ke,
                ve / q.locationSize * we * ke,
                Ye
              );
          }
        } else if (ee !== void 0) {
          const me = ee[I];
          if (me !== void 0)
            switch (me.length) {
              case 2:
                i.vertexAttrib2fv(q.location, me);
                break;
              case 3:
                i.vertexAttrib3fv(q.location, me);
                break;
              case 4:
                i.vertexAttrib4fv(q.location, me);
                break;
              default:
                i.vertexAttrib1fv(q.location, me);
            }
        }
      }
    }
    L();
  }
  function M() {
    ce();
    for (const H in a) {
      const G = a[H];
      for (const Q in G) {
        const X = G[Q];
        for (const j in X)
          g(X[j].object), delete X[j];
        delete G[Q];
      }
      delete a[H];
    }
  }
  function y(H) {
    if (a[H.id] === void 0)
      return;
    const G = a[H.id];
    for (const Q in G) {
      const X = G[Q];
      for (const j in X)
        g(X[j].object), delete X[j];
      delete G[Q];
    }
    delete a[H.id];
  }
  function Y(H) {
    for (const G in a) {
      const Q = a[G];
      if (Q[H.id] === void 0)
        continue;
      const X = Q[H.id];
      for (const j in X)
        g(X[j].object), delete X[j];
      delete Q[H.id];
    }
  }
  function ce() {
    B(), h = !0, c !== l && (c = l, m(c.object));
  }
  function B() {
    l.geometry = null, l.program = null, l.wireframe = !1;
  }
  return {
    setup: f,
    reset: ce,
    resetDefaultState: B,
    dispose: M,
    releaseStatesOfGeometry: y,
    releaseStatesOfProgram: Y,
    initAttributes: _,
    enableAttribute: T,
    disableUnusedAttributes: L
  };
}
function gf(i, e, t, n) {
  const r = n.isWebGL2;
  let s;
  function o(c) {
    s = c;
  }
  function a(c, h) {
    i.drawArrays(s, c, h), t.update(h, s, 1);
  }
  function l(c, h, f) {
    if (f === 0)
      return;
    let u, m;
    if (r)
      u = i, m = "drawArraysInstanced";
    else if (u = e.get("ANGLE_instanced_arrays"), m = "drawArraysInstancedANGLE", u === null) {
      console.error("THREE.WebGLBufferRenderer: using THREE.InstancedBufferGeometry but hardware does not support extension ANGLE_instanced_arrays.");
      return;
    }
    u[m](s, c, h, f), t.update(h, s, f);
  }
  this.setMode = o, this.render = a, this.renderInstances = l;
}
function _f(i, e, t) {
  let n;
  function r() {
    if (n !== void 0)
      return n;
    if (e.has("EXT_texture_filter_anisotropic") === !0) {
      const w = e.get("EXT_texture_filter_anisotropic");
      n = i.getParameter(w.MAX_TEXTURE_MAX_ANISOTROPY_EXT);
    } else
      n = 0;
    return n;
  }
  function s(w) {
    if (w === "highp") {
      if (i.getShaderPrecisionFormat(i.VERTEX_SHADER, i.HIGH_FLOAT).precision > 0 && i.getShaderPrecisionFormat(i.FRAGMENT_SHADER, i.HIGH_FLOAT).precision > 0)
        return "highp";
      w = "mediump";
    }
    return w === "mediump" && i.getShaderPrecisionFormat(i.VERTEX_SHADER, i.MEDIUM_FLOAT).precision > 0 && i.getShaderPrecisionFormat(i.FRAGMENT_SHADER, i.MEDIUM_FLOAT).precision > 0 ? "mediump" : "lowp";
  }
  const o = typeof WebGL2RenderingContext < "u" && i.constructor.name === "WebGL2RenderingContext";
  let a = t.precision !== void 0 ? t.precision : "highp";
  const l = s(a);
  l !== a && (console.warn("THREE.WebGLRenderer:", a, "not supported, using", l, "instead."), a = l);
  const c = o || e.has("WEBGL_draw_buffers"), h = t.logarithmicDepthBuffer === !0, f = i.getParameter(i.MAX_TEXTURE_IMAGE_UNITS), u = i.getParameter(i.MAX_VERTEX_TEXTURE_IMAGE_UNITS), m = i.getParameter(i.MAX_TEXTURE_SIZE), g = i.getParameter(i.MAX_CUBE_MAP_TEXTURE_SIZE), x = i.getParameter(i.MAX_VERTEX_ATTRIBS), p = i.getParameter(i.MAX_VERTEX_UNIFORM_VECTORS), d = i.getParameter(i.MAX_VARYING_VECTORS), A = i.getParameter(i.MAX_FRAGMENT_UNIFORM_VECTORS), _ = u > 0, T = o || e.has("OES_texture_float"), C = _ && T, L = o ? i.getParameter(i.MAX_SAMPLES) : 0;
  return {
    isWebGL2: o,
    drawBuffers: c,
    getMaxAnisotropy: r,
    getMaxPrecision: s,
    precision: a,
    logarithmicDepthBuffer: h,
    maxTextures: f,
    maxVertexTextures: u,
    maxTextureSize: m,
    maxCubemapSize: g,
    maxAttributes: x,
    maxVertexUniforms: p,
    maxVaryings: d,
    maxFragmentUniforms: A,
    vertexTextures: _,
    floatFragmentTextures: T,
    floatVertexTextures: C,
    maxSamples: L
  };
}
function vf(i) {
  const e = this;
  let t = null, n = 0, r = !1, s = !1;
  const o = new en(), a = new Be(), l = { value: null, needsUpdate: !1 };
  this.uniform = l, this.numPlanes = 0, this.numIntersection = 0, this.init = function(f, u) {
    const m = f.length !== 0 || u || // enable state of previous frame - the clipping code has to
    // run another frame in order to reset the state:
    n !== 0 || r;
    return r = u, n = f.length, m;
  }, this.beginShadows = function() {
    s = !0, h(null);
  }, this.endShadows = function() {
    s = !1;
  }, this.setGlobalState = function(f, u) {
    t = h(f, u, 0);
  }, this.setState = function(f, u, m) {
    const g = f.clippingPlanes, x = f.clipIntersection, p = f.clipShadows, d = i.get(f);
    if (!r || g === null || g.length === 0 || s && !p)
      s ? h(null) : c();
    else {
      const A = s ? 0 : n, _ = A * 4;
      let T = d.clippingState || null;
      l.value = T, T = h(g, u, _, m);
      for (let C = 0; C !== _; ++C)
        T[C] = t[C];
      d.clippingState = T, this.numIntersection = x ? this.numPlanes : 0, this.numPlanes += A;
    }
  };
  function c() {
    l.value !== t && (l.value = t, l.needsUpdate = n > 0), e.numPlanes = n, e.numIntersection = 0;
  }
  function h(f, u, m, g) {
    const x = f !== null ? f.length : 0;
    let p = null;
    if (x !== 0) {
      if (p = l.value, g !== !0 || p === null) {
        const d = m + x * 4, A = u.matrixWorldInverse;
        a.getNormalMatrix(A), (p === null || p.length < d) && (p = new Float32Array(d));
        for (let _ = 0, T = m; _ !== x; ++_, T += 4)
          o.copy(f[_]).applyMatrix4(A, a), o.normal.toArray(p, T), p[T + 3] = o.constant;
      }
      l.value = p, l.needsUpdate = !0;
    }
    return e.numPlanes = x, e.numIntersection = 0, p;
  }
}
function xf(i) {
  let e = /* @__PURE__ */ new WeakMap();
  function t(o, a) {
    return a === qr ? o.mapping = ei : a === jr && (o.mapping = ti), o;
  }
  function n(o) {
    if (o && o.isTexture && o.isRenderTargetTexture === !1) {
      const a = o.mapping;
      if (a === qr || a === jr)
        if (e.has(o)) {
          const l = e.get(o).texture;
          return t(l, o.mapping);
        } else {
          const l = o.image;
          if (l && l.height > 0) {
            const c = new Uc(l.height / 2);
            return c.fromEquirectangularTexture(i, o), e.set(o, c), o.addEventListener("dispose", r), t(c.texture, o.mapping);
          } else
            return null;
        }
    }
    return o;
  }
  function r(o) {
    const a = o.target;
    a.removeEventListener("dispose", r);
    const l = e.get(a);
    l !== void 0 && (e.delete(a), l.dispose());
  }
  function s() {
    e = /* @__PURE__ */ new WeakMap();
  }
  return {
    get: n,
    dispose: s
  };
}
class Ao extends Eo {
  constructor(e = -1, t = 1, n = 1, r = -1, s = 0.1, o = 2e3) {
    super(), this.isOrthographicCamera = !0, this.type = "OrthographicCamera", this.zoom = 1, this.view = null, this.left = e, this.right = t, this.top = n, this.bottom = r, this.near = s, this.far = o, this.updateProjectionMatrix();
  }
  copy(e, t) {
    return super.copy(e, t), this.left = e.left, this.right = e.right, this.top = e.top, this.bottom = e.bottom, this.near = e.near, this.far = e.far, this.zoom = e.zoom, this.view = e.view === null ? null : Object.assign({}, e.view), this;
  }
  setViewOffset(e, t, n, r, s, o) {
    this.view === null && (this.view = {
      enabled: !0,
      fullWidth: 1,
      fullHeight: 1,
      offsetX: 0,
      offsetY: 0,
      width: 1,
      height: 1
    }), this.view.enabled = !0, this.view.fullWidth = e, this.view.fullHeight = t, this.view.offsetX = n, this.view.offsetY = r, this.view.width = s, this.view.height = o, this.updateProjectionMatrix();
  }
  clearViewOffset() {
    this.view !== null && (this.view.enabled = !1), this.updateProjectionMatrix();
  }
  updateProjectionMatrix() {
    const e = (this.right - this.left) / (2 * this.zoom), t = (this.top - this.bottom) / (2 * this.zoom), n = (this.right + this.left) / 2, r = (this.top + this.bottom) / 2;
    let s = n - e, o = n + e, a = r + t, l = r - t;
    if (this.view !== null && this.view.enabled) {
      const c = (this.right - this.left) / this.view.fullWidth / this.zoom, h = (this.top - this.bottom) / this.view.fullHeight / this.zoom;
      s += c * this.view.offsetX, o = s + c * this.view.width, a -= h * this.view.offsetY, l = a - h * this.view.height;
    }
    this.projectionMatrix.makeOrthographic(s, o, a, l, this.near, this.far, this.coordinateSystem), this.projectionMatrixInverse.copy(this.projectionMatrix).invert();
  }
  toJSON(e) {
    const t = super.toJSON(e);
    return t.object.zoom = this.zoom, t.object.left = this.left, t.object.right = this.right, t.object.top = this.top, t.object.bottom = this.bottom, t.object.near = this.near, t.object.far = this.far, this.view !== null && (t.object.view = Object.assign({}, this.view)), t;
  }
}
const Kn = 4, _a = [0.125, 0.215, 0.35, 0.446, 0.526, 0.582], vn = 20, Fr = /* @__PURE__ */ new Ao(), va = /* @__PURE__ */ new We();
let Br = null;
const gn = (1 + Math.sqrt(5)) / 2, qn = 1 / gn, xa = [
  /* @__PURE__ */ new U(1, 1, 1),
  /* @__PURE__ */ new U(-1, 1, 1),
  /* @__PURE__ */ new U(1, 1, -1),
  /* @__PURE__ */ new U(-1, 1, -1),
  /* @__PURE__ */ new U(0, gn, qn),
  /* @__PURE__ */ new U(0, gn, -qn),
  /* @__PURE__ */ new U(qn, 0, gn),
  /* @__PURE__ */ new U(-qn, 0, gn),
  /* @__PURE__ */ new U(gn, qn, 0),
  /* @__PURE__ */ new U(-gn, qn, 0)
];
class Ma {
  constructor(e) {
    this._renderer = e, this._pingPongRenderTarget = null, this._lodMax = 0, this._cubeSize = 0, this._lodPlanes = [], this._sizeLods = [], this._sigmas = [], this._blurMaterial = null, this._cubemapMaterial = null, this._equirectMaterial = null, this._compileMaterial(this._blurMaterial);
  }
  /**
   * Generates a PMREM from a supplied Scene, which can be faster than using an
   * image if networking bandwidth is low. Optional sigma specifies a blur radius
   * in radians to be applied to the scene before PMREM generation. Optional near
   * and far planes ensure the scene is rendered in its entirety (the cubeCamera
   * is placed at the origin).
   */
  fromScene(e, t = 0, n = 0.1, r = 100) {
    Br = this._renderer.getRenderTarget(), this._setSize(256);
    const s = this._allocateTargets();
    return s.depthBuffer = !0, this._sceneToCubeUV(e, n, r, s), t > 0 && this._blur(s, 0, 0, t), this._applyPMREM(s), this._cleanup(s), s;
  }
  /**
   * Generates a PMREM from an equirectangular texture, which can be either LDR
   * or HDR. The ideal input image size is 1k (1024 x 512),
   * as this matches best with the 256 x 256 cubemap output.
   */
  fromEquirectangular(e, t = null) {
    return this._fromTexture(e, t);
  }
  /**
   * Generates a PMREM from an cubemap texture, which can be either LDR
   * or HDR. The ideal input cube size is 256 x 256,
   * as this matches best with the 256 x 256 cubemap output.
   */
  fromCubemap(e, t = null) {
    return this._fromTexture(e, t);
  }
  /**
   * Pre-compiles the cubemap shader. You can get faster start-up by invoking this method during
   * your texture's network fetch for increased concurrency.
   */
  compileCubemapShader() {
    this._cubemapMaterial === null && (this._cubemapMaterial = ya(), this._compileMaterial(this._cubemapMaterial));
  }
  /**
   * Pre-compiles the equirectangular shader. You can get faster start-up by invoking this method during
   * your texture's network fetch for increased concurrency.
   */
  compileEquirectangularShader() {
    this._equirectMaterial === null && (this._equirectMaterial = Ea(), this._compileMaterial(this._equirectMaterial));
  }
  /**
   * Disposes of the PMREMGenerator's internal memory. Note that PMREMGenerator is a static class,
   * so you should not need more than one PMREMGenerator object. If you do, calling dispose() on
   * one of them will cause any others to also become unusable.
   */
  dispose() {
    this._dispose(), this._cubemapMaterial !== null && this._cubemapMaterial.dispose(), this._equirectMaterial !== null && this._equirectMaterial.dispose();
  }
  // private interface
  _setSize(e) {
    this._lodMax = Math.floor(Math.log2(e)), this._cubeSize = Math.pow(2, this._lodMax);
  }
  _dispose() {
    this._blurMaterial !== null && this._blurMaterial.dispose(), this._pingPongRenderTarget !== null && this._pingPongRenderTarget.dispose();
    for (let e = 0; e < this._lodPlanes.length; e++)
      this._lodPlanes[e].dispose();
  }
  _cleanup(e) {
    this._renderer.setRenderTarget(Br), e.scissorTest = !1, ji(e, 0, 0, e.width, e.height);
  }
  _fromTexture(e, t) {
    e.mapping === ei || e.mapping === ti ? this._setSize(e.image.length === 0 ? 16 : e.image[0].width || e.image[0].image.width) : this._setSize(e.image.width / 4), Br = this._renderer.getRenderTarget();
    const n = t || this._allocateTargets();
    return this._textureToCubeUV(e, n), this._applyPMREM(n), this._cleanup(n), n;
  }
  _allocateTargets() {
    const e = 3 * Math.max(this._cubeSize, 112), t = 4 * this._cubeSize, n = {
      magFilter: Tt,
      minFilter: Tt,
      generateMipmaps: !1,
      type: xi,
      format: Ut,
      colorSpace: Ft,
      depthBuffer: !1
    }, r = Sa(e, t, n);
    if (this._pingPongRenderTarget === null || this._pingPongRenderTarget.width !== e || this._pingPongRenderTarget.height !== t) {
      this._pingPongRenderTarget !== null && this._dispose(), this._pingPongRenderTarget = Sa(e, t, n);
      const { _lodMax: s } = this;
      ({ sizeLods: this._sizeLods, lodPlanes: this._lodPlanes, sigmas: this._sigmas } = Mf(s)), this._blurMaterial = Sf(s, e, t);
    }
    return r;
  }
  _compileMaterial(e) {
    const t = new Nt(this._lodPlanes[0], e);
    this._renderer.compile(t, Fr);
  }
  _sceneToCubeUV(e, t, n, r) {
    const a = new At(90, 1, t, n), l = [1, -1, 1, 1, 1, 1], c = [1, 1, 1, -1, -1, -1], h = this._renderer, f = h.autoClear, u = h.toneMapping;
    h.getClearColor(va), h.toneMapping = an, h.autoClear = !1;
    const m = new vo({
      name: "PMREM.Background",
      side: gt,
      depthWrite: !1,
      depthTest: !1
    }), g = new Nt(new bi(), m);
    let x = !1;
    const p = e.background;
    p ? p.isColor && (m.color.copy(p), e.background = null, x = !0) : (m.color.copy(va), x = !0);
    for (let d = 0; d < 6; d++) {
      const A = d % 3;
      A === 0 ? (a.up.set(0, l[d], 0), a.lookAt(c[d], 0, 0)) : A === 1 ? (a.up.set(0, 0, l[d]), a.lookAt(0, c[d], 0)) : (a.up.set(0, l[d], 0), a.lookAt(0, 0, c[d]));
      const _ = this._cubeSize;
      ji(r, A * _, d > 2 ? _ : 0, _, _), h.setRenderTarget(r), x && h.render(g, a), h.render(e, a);
    }
    g.geometry.dispose(), g.material.dispose(), h.toneMapping = u, h.autoClear = f, e.background = p;
  }
  _textureToCubeUV(e, t) {
    const n = this._renderer, r = e.mapping === ei || e.mapping === ti;
    r ? (this._cubemapMaterial === null && (this._cubemapMaterial = ya()), this._cubemapMaterial.uniforms.flipEnvMap.value = e.isRenderTargetTexture === !1 ? -1 : 1) : this._equirectMaterial === null && (this._equirectMaterial = Ea());
    const s = r ? this._cubemapMaterial : this._equirectMaterial, o = new Nt(this._lodPlanes[0], s), a = s.uniforms;
    a.envMap.value = e;
    const l = this._cubeSize;
    ji(t, 0, 0, 3 * l, 2 * l), n.setRenderTarget(t), n.render(o, Fr);
  }
  _applyPMREM(e) {
    const t = this._renderer, n = t.autoClear;
    t.autoClear = !1;
    for (let r = 1; r < this._lodPlanes.length; r++) {
      const s = Math.sqrt(this._sigmas[r] * this._sigmas[r] - this._sigmas[r - 1] * this._sigmas[r - 1]), o = xa[(r - 1) % xa.length];
      this._blur(e, r - 1, r, s, o);
    }
    t.autoClear = n;
  }
  /**
   * This is a two-pass Gaussian blur for a cubemap. Normally this is done
   * vertically and horizontally, but this breaks down on a cube. Here we apply
   * the blur latitudinally (around the poles), and then longitudinally (towards
   * the poles) to approximate the orthogonally-separable blur. It is least
   * accurate at the poles, but still does a decent job.
   */
  _blur(e, t, n, r, s) {
    const o = this._pingPongRenderTarget;
    this._halfBlur(
      e,
      o,
      t,
      n,
      r,
      "latitudinal",
      s
    ), this._halfBlur(
      o,
      e,
      n,
      n,
      r,
      "longitudinal",
      s
    );
  }
  _halfBlur(e, t, n, r, s, o, a) {
    const l = this._renderer, c = this._blurMaterial;
    o !== "latitudinal" && o !== "longitudinal" && console.error(
      "blur direction must be either latitudinal or longitudinal!"
    );
    const h = 3, f = new Nt(this._lodPlanes[r], c), u = c.uniforms, m = this._sizeLods[n] - 1, g = isFinite(s) ? Math.PI / (2 * m) : 2 * Math.PI / (2 * vn - 1), x = s / g, p = isFinite(s) ? 1 + Math.floor(h * x) : vn;
    p > vn && console.warn(`sigmaRadians, ${s}, is too large and will clip, as it requested ${p} samples when the maximum is set to ${vn}`);
    const d = [];
    let A = 0;
    for (let w = 0; w < vn; ++w) {
      const V = w / x, M = Math.exp(-V * V / 2);
      d.push(M), w === 0 ? A += M : w < p && (A += 2 * M);
    }
    for (let w = 0; w < d.length; w++)
      d[w] = d[w] / A;
    u.envMap.value = e.texture, u.samples.value = p, u.weights.value = d, u.latitudinal.value = o === "latitudinal", a && (u.poleAxis.value = a);
    const { _lodMax: _ } = this;
    u.dTheta.value = g, u.mipInt.value = _ - n;
    const T = this._sizeLods[r], C = 3 * T * (r > _ - Kn ? r - _ + Kn : 0), L = 4 * (this._cubeSize - T);
    ji(t, C, L, 3 * T, 2 * T), l.setRenderTarget(t), l.render(f, Fr);
  }
}
function Mf(i) {
  const e = [], t = [], n = [];
  let r = i;
  const s = i - Kn + 1 + _a.length;
  for (let o = 0; o < s; o++) {
    const a = Math.pow(2, r);
    t.push(a);
    let l = 1 / a;
    o > i - Kn ? l = _a[o - i + Kn - 1] : o === 0 && (l = 0), n.push(l);
    const c = 1 / (a - 2), h = -c, f = 1 + c, u = [h, h, f, h, f, f, h, h, f, f, h, f], m = 6, g = 6, x = 3, p = 2, d = 1, A = new Float32Array(x * g * m), _ = new Float32Array(p * g * m), T = new Float32Array(d * g * m);
    for (let L = 0; L < m; L++) {
      const w = L % 3 * 2 / 3 - 1, V = L > 2 ? 0 : -1, M = [
        w,
        V,
        0,
        w + 2 / 3,
        V,
        0,
        w + 2 / 3,
        V + 1,
        0,
        w,
        V,
        0,
        w + 2 / 3,
        V + 1,
        0,
        w,
        V + 1,
        0
      ];
      A.set(M, x * g * L), _.set(u, p * g * L);
      const y = [L, L, L, L, L, L];
      T.set(y, d * g * L);
    }
    const C = new cn();
    C.setAttribute("position", new Ot(A, x)), C.setAttribute("uv", new Ot(_, p)), C.setAttribute("faceIndex", new Ot(T, d)), e.push(C), r > Kn && r--;
  }
  return { lodPlanes: e, sizeLods: t, sigmas: n };
}
function Sa(i, e, t) {
  const n = new Tn(i, e, t);
  return n.texture.mapping = rr, n.texture.name = "PMREM.cubeUv", n.scissorTest = !0, n;
}
function ji(i, e, t, n, r) {
  i.viewport.set(e, t, n, r), i.scissor.set(e, t, n, r);
}
function Sf(i, e, t) {
  const n = new Float32Array(vn), r = new U(0, 1, 0);
  return new bn({
    name: "SphericalGaussianBlur",
    defines: {
      n: vn,
      CUBEUV_TEXEL_WIDTH: 1 / e,
      CUBEUV_TEXEL_HEIGHT: 1 / t,
      CUBEUV_MAX_MIP: `${i}.0`
    },
    uniforms: {
      envMap: { value: null },
      samples: { value: 1 },
      weights: { value: n },
      latitudinal: { value: !1 },
      dTheta: { value: 0 },
      mipInt: { value: 0 },
      poleAxis: { value: r }
    },
    vertexShader: fs(),
    fragmentShader: (
      /* glsl */
      `

			precision mediump float;
			precision mediump int;

			varying vec3 vOutputDirection;

			uniform sampler2D envMap;
			uniform int samples;
			uniform float weights[ n ];
			uniform bool latitudinal;
			uniform float dTheta;
			uniform float mipInt;
			uniform vec3 poleAxis;

			#define ENVMAP_TYPE_CUBE_UV
			#include <cube_uv_reflection_fragment>

			vec3 getSample( float theta, vec3 axis ) {

				float cosTheta = cos( theta );
				// Rodrigues' axis-angle rotation
				vec3 sampleDirection = vOutputDirection * cosTheta
					+ cross( axis, vOutputDirection ) * sin( theta )
					+ axis * dot( axis, vOutputDirection ) * ( 1.0 - cosTheta );

				return bilinearCubeUV( envMap, sampleDirection, mipInt );

			}

			void main() {

				vec3 axis = latitudinal ? poleAxis : cross( poleAxis, vOutputDirection );

				if ( all( equal( axis, vec3( 0.0 ) ) ) ) {

					axis = vec3( vOutputDirection.z, 0.0, - vOutputDirection.x );

				}

				axis = normalize( axis );

				gl_FragColor = vec4( 0.0, 0.0, 0.0, 1.0 );
				gl_FragColor.rgb += weights[ 0 ] * getSample( 0.0, axis );

				for ( int i = 1; i < n; i++ ) {

					if ( i >= samples ) {

						break;

					}

					float theta = dTheta * float( i );
					gl_FragColor.rgb += weights[ i ] * getSample( -1.0 * theta, axis );
					gl_FragColor.rgb += weights[ i ] * getSample( theta, axis );

				}

			}
		`
    ),
    blending: sn,
    depthTest: !1,
    depthWrite: !1
  });
}
function Ea() {
  return new bn({
    name: "EquirectangularToCubeUV",
    uniforms: {
      envMap: { value: null }
    },
    vertexShader: fs(),
    fragmentShader: (
      /* glsl */
      `

			precision mediump float;
			precision mediump int;

			varying vec3 vOutputDirection;

			uniform sampler2D envMap;

			#include <common>

			void main() {

				vec3 outputDirection = normalize( vOutputDirection );
				vec2 uv = equirectUv( outputDirection );

				gl_FragColor = vec4( texture2D ( envMap, uv ).rgb, 1.0 );

			}
		`
    ),
    blending: sn,
    depthTest: !1,
    depthWrite: !1
  });
}
function ya() {
  return new bn({
    name: "CubemapToCubeUV",
    uniforms: {
      envMap: { value: null },
      flipEnvMap: { value: -1 }
    },
    vertexShader: fs(),
    fragmentShader: (
      /* glsl */
      `

			precision mediump float;
			precision mediump int;

			uniform float flipEnvMap;

			varying vec3 vOutputDirection;

			uniform samplerCube envMap;

			void main() {

				gl_FragColor = textureCube( envMap, vec3( flipEnvMap * vOutputDirection.x, vOutputDirection.yz ) );

			}
		`
    ),
    blending: sn,
    depthTest: !1,
    depthWrite: !1
  });
}
function fs() {
  return (
    /* glsl */
    `

		precision mediump float;
		precision mediump int;

		attribute float faceIndex;

		varying vec3 vOutputDirection;

		// RH coordinate system; PMREM face-indexing convention
		vec3 getDirection( vec2 uv, float face ) {

			uv = 2.0 * uv - 1.0;

			vec3 direction = vec3( uv, 1.0 );

			if ( face == 0.0 ) {

				direction = direction.zyx; // ( 1, v, u ) pos x

			} else if ( face == 1.0 ) {

				direction = direction.xzy;
				direction.xz *= -1.0; // ( -u, 1, -v ) pos y

			} else if ( face == 2.0 ) {

				direction.x *= -1.0; // ( -u, v, 1 ) pos z

			} else if ( face == 3.0 ) {

				direction = direction.zyx;
				direction.xz *= -1.0; // ( -1, v, -u ) neg x

			} else if ( face == 4.0 ) {

				direction = direction.xzy;
				direction.xy *= -1.0; // ( -u, -1, v ) neg y

			} else if ( face == 5.0 ) {

				direction.z *= -1.0; // ( u, v, -1 ) neg z

			}

			return direction;

		}

		void main() {

			vOutputDirection = getDirection( uv, faceIndex );
			gl_Position = vec4( position, 1.0 );

		}
	`
  );
}
function Ef(i) {
  let e = /* @__PURE__ */ new WeakMap(), t = null;
  function n(a) {
    if (a && a.isTexture) {
      const l = a.mapping, c = l === qr || l === jr, h = l === ei || l === ti;
      if (c || h)
        if (a.isRenderTargetTexture && a.needsPMREMUpdate === !0) {
          a.needsPMREMUpdate = !1;
          let f = e.get(a);
          return t === null && (t = new Ma(i)), f = c ? t.fromEquirectangular(a, f) : t.fromCubemap(a, f), e.set(a, f), f.texture;
        } else {
          if (e.has(a))
            return e.get(a).texture;
          {
            const f = a.image;
            if (c && f && f.height > 0 || h && f && r(f)) {
              t === null && (t = new Ma(i));
              const u = c ? t.fromEquirectangular(a) : t.fromCubemap(a);
              return e.set(a, u), a.addEventListener("dispose", s), u.texture;
            } else
              return null;
          }
        }
    }
    return a;
  }
  function r(a) {
    let l = 0;
    const c = 6;
    for (let h = 0; h < c; h++)
      a[h] !== void 0 && l++;
    return l === c;
  }
  function s(a) {
    const l = a.target;
    l.removeEventListener("dispose", s);
    const c = e.get(l);
    c !== void 0 && (e.delete(l), c.dispose());
  }
  function o() {
    e = /* @__PURE__ */ new WeakMap(), t !== null && (t.dispose(), t = null);
  }
  return {
    get: n,
    dispose: o
  };
}
function yf(i) {
  const e = {};
  function t(n) {
    if (e[n] !== void 0)
      return e[n];
    let r;
    switch (n) {
      case "WEBGL_depth_texture":
        r = i.getExtension("WEBGL_depth_texture") || i.getExtension("MOZ_WEBGL_depth_texture") || i.getExtension("WEBKIT_WEBGL_depth_texture");
        break;
      case "EXT_texture_filter_anisotropic":
        r = i.getExtension("EXT_texture_filter_anisotropic") || i.getExtension("MOZ_EXT_texture_filter_anisotropic") || i.getExtension("WEBKIT_EXT_texture_filter_anisotropic");
        break;
      case "WEBGL_compressed_texture_s3tc":
        r = i.getExtension("WEBGL_compressed_texture_s3tc") || i.getExtension("MOZ_WEBGL_compressed_texture_s3tc") || i.getExtension("WEBKIT_WEBGL_compressed_texture_s3tc");
        break;
      case "WEBGL_compressed_texture_pvrtc":
        r = i.getExtension("WEBGL_compressed_texture_pvrtc") || i.getExtension("WEBKIT_WEBGL_compressed_texture_pvrtc");
        break;
      default:
        r = i.getExtension(n);
    }
    return e[n] = r, r;
  }
  return {
    has: function(n) {
      return t(n) !== null;
    },
    init: function(n) {
      n.isWebGL2 ? t("EXT_color_buffer_float") : (t("WEBGL_depth_texture"), t("OES_texture_float"), t("OES_texture_half_float"), t("OES_texture_half_float_linear"), t("OES_standard_derivatives"), t("OES_element_index_uint"), t("OES_vertex_array_object"), t("ANGLE_instanced_arrays")), t("OES_texture_float_linear"), t("EXT_color_buffer_half_float"), t("WEBGL_multisampled_render_to_texture");
    },
    get: function(n) {
      const r = t(n);
      return r === null && console.warn("THREE.WebGLRenderer: " + n + " extension not supported."), r;
    }
  };
}
function Tf(i, e, t, n) {
  const r = {}, s = /* @__PURE__ */ new WeakMap();
  function o(f) {
    const u = f.target;
    u.index !== null && e.remove(u.index);
    for (const g in u.attributes)
      e.remove(u.attributes[g]);
    for (const g in u.morphAttributes) {
      const x = u.morphAttributes[g];
      for (let p = 0, d = x.length; p < d; p++)
        e.remove(x[p]);
    }
    u.removeEventListener("dispose", o), delete r[u.id];
    const m = s.get(u);
    m && (e.remove(m), s.delete(u)), n.releaseStatesOfGeometry(u), u.isInstancedBufferGeometry === !0 && delete u._maxInstanceCount, t.memory.geometries--;
  }
  function a(f, u) {
    return r[u.id] === !0 || (u.addEventListener("dispose", o), r[u.id] = !0, t.memory.geometries++), u;
  }
  function l(f) {
    const u = f.attributes;
    for (const g in u)
      e.update(u[g], i.ARRAY_BUFFER);
    const m = f.morphAttributes;
    for (const g in m) {
      const x = m[g];
      for (let p = 0, d = x.length; p < d; p++)
        e.update(x[p], i.ARRAY_BUFFER);
    }
  }
  function c(f) {
    const u = [], m = f.index, g = f.attributes.position;
    let x = 0;
    if (m !== null) {
      const A = m.array;
      x = m.version;
      for (let _ = 0, T = A.length; _ < T; _ += 3) {
        const C = A[_ + 0], L = A[_ + 1], w = A[_ + 2];
        u.push(C, L, L, w, w, C);
      }
    } else if (g !== void 0) {
      const A = g.array;
      x = g.version;
      for (let _ = 0, T = A.length / 3 - 1; _ < T; _ += 3) {
        const C = _ + 0, L = _ + 1, w = _ + 2;
        u.push(C, L, L, w, w, C);
      }
    } else
      return;
    const p = new (fo(u) ? Mo : xo)(u, 1);
    p.version = x;
    const d = s.get(f);
    d && e.remove(d), s.set(f, p);
  }
  function h(f) {
    const u = s.get(f);
    if (u) {
      const m = f.index;
      m !== null && u.version < m.version && c(f);
    } else
      c(f);
    return s.get(f);
  }
  return {
    get: a,
    update: l,
    getWireframeAttribute: h
  };
}
function Af(i, e, t, n) {
  const r = n.isWebGL2;
  let s;
  function o(u) {
    s = u;
  }
  let a, l;
  function c(u) {
    a = u.type, l = u.bytesPerElement;
  }
  function h(u, m) {
    i.drawElements(s, m, a, u * l), t.update(m, s, 1);
  }
  function f(u, m, g) {
    if (g === 0)
      return;
    let x, p;
    if (r)
      x = i, p = "drawElementsInstanced";
    else if (x = e.get("ANGLE_instanced_arrays"), p = "drawElementsInstancedANGLE", x === null) {
      console.error("THREE.WebGLIndexedBufferRenderer: using THREE.InstancedBufferGeometry but hardware does not support extension ANGLE_instanced_arrays.");
      return;
    }
    x[p](s, m, a, u * l, g), t.update(m, s, g);
  }
  this.setMode = o, this.setIndex = c, this.render = h, this.renderInstances = f;
}
function bf(i) {
  const e = {
    geometries: 0,
    textures: 0
  }, t = {
    frame: 0,
    calls: 0,
    triangles: 0,
    points: 0,
    lines: 0
  };
  function n(s, o, a) {
    switch (t.calls++, o) {
      case i.TRIANGLES:
        t.triangles += a * (s / 3);
        break;
      case i.LINES:
        t.lines += a * (s / 2);
        break;
      case i.LINE_STRIP:
        t.lines += a * (s - 1);
        break;
      case i.LINE_LOOP:
        t.lines += a * s;
        break;
      case i.POINTS:
        t.points += a * s;
        break;
      default:
        console.error("THREE.WebGLInfo: Unknown draw mode:", o);
        break;
    }
  }
  function r() {
    t.calls = 0, t.triangles = 0, t.points = 0, t.lines = 0;
  }
  return {
    memory: e,
    render: t,
    programs: null,
    autoReset: !0,
    reset: r,
    update: n
  };
}
function wf(i, e) {
  return i[0] - e[0];
}
function Rf(i, e) {
  return Math.abs(e[1]) - Math.abs(i[1]);
}
function Cf(i, e, t) {
  const n = {}, r = new Float32Array(8), s = /* @__PURE__ */ new WeakMap(), o = new ot(), a = [];
  for (let c = 0; c < 8; c++)
    a[c] = [c, 0];
  function l(c, h, f) {
    const u = c.morphTargetInfluences;
    if (e.isWebGL2 === !0) {
      const m = h.morphAttributes.position || h.morphAttributes.normal || h.morphAttributes.color, g = m !== void 0 ? m.length : 0;
      let x = s.get(h);
      if (x === void 0 || x.count !== g) {
        let H = function() {
          ce.dispose(), s.delete(h), h.removeEventListener("dispose", H);
        };
        x !== void 0 && x.texture.dispose();
        const A = h.morphAttributes.position !== void 0, _ = h.morphAttributes.normal !== void 0, T = h.morphAttributes.color !== void 0, C = h.morphAttributes.position || [], L = h.morphAttributes.normal || [], w = h.morphAttributes.color || [];
        let V = 0;
        A === !0 && (V = 1), _ === !0 && (V = 2), T === !0 && (V = 3);
        let M = h.attributes.position.count * V, y = 1;
        M > e.maxTextureSize && (y = Math.ceil(M / e.maxTextureSize), M = e.maxTextureSize);
        const Y = new Float32Array(M * y * 4 * g), ce = new go(Y, M, y, g);
        ce.type = rn, ce.needsUpdate = !0;
        const B = V * 4;
        for (let G = 0; G < g; G++) {
          const Q = C[G], X = L[G], j = w[G], J = M * y * 4 * G;
          for (let ee = 0; ee < Q.count; ee++) {
            const I = ee * B;
            A === !0 && (o.fromBufferAttribute(Q, ee), Y[J + I + 0] = o.x, Y[J + I + 1] = o.y, Y[J + I + 2] = o.z, Y[J + I + 3] = 0), _ === !0 && (o.fromBufferAttribute(X, ee), Y[J + I + 4] = o.x, Y[J + I + 5] = o.y, Y[J + I + 6] = o.z, Y[J + I + 7] = 0), T === !0 && (o.fromBufferAttribute(j, ee), Y[J + I + 8] = o.x, Y[J + I + 9] = o.y, Y[J + I + 10] = o.z, Y[J + I + 11] = j.itemSize === 4 ? o.w : 1);
          }
        }
        x = {
          count: g,
          texture: ce,
          size: new oe(M, y)
        }, s.set(h, x), h.addEventListener("dispose", H);
      }
      let p = 0;
      for (let A = 0; A < u.length; A++)
        p += u[A];
      const d = h.morphTargetsRelative ? 1 : 1 - p;
      f.getUniforms().setValue(i, "morphTargetBaseInfluence", d), f.getUniforms().setValue(i, "morphTargetInfluences", u), f.getUniforms().setValue(i, "morphTargetsTexture", x.texture, t), f.getUniforms().setValue(i, "morphTargetsTextureSize", x.size);
    } else {
      const m = u === void 0 ? 0 : u.length;
      let g = n[h.id];
      if (g === void 0 || g.length !== m) {
        g = [];
        for (let _ = 0; _ < m; _++)
          g[_] = [_, 0];
        n[h.id] = g;
      }
      for (let _ = 0; _ < m; _++) {
        const T = g[_];
        T[0] = _, T[1] = u[_];
      }
      g.sort(Rf);
      for (let _ = 0; _ < 8; _++)
        _ < m && g[_][1] ? (a[_][0] = g[_][0], a[_][1] = g[_][1]) : (a[_][0] = Number.MAX_SAFE_INTEGER, a[_][1] = 0);
      a.sort(wf);
      const x = h.morphAttributes.position, p = h.morphAttributes.normal;
      let d = 0;
      for (let _ = 0; _ < 8; _++) {
        const T = a[_], C = T[0], L = T[1];
        C !== Number.MAX_SAFE_INTEGER && L ? (x && h.getAttribute("morphTarget" + _) !== x[C] && h.setAttribute("morphTarget" + _, x[C]), p && h.getAttribute("morphNormal" + _) !== p[C] && h.setAttribute("morphNormal" + _, p[C]), r[_] = L, d += L) : (x && h.hasAttribute("morphTarget" + _) === !0 && h.deleteAttribute("morphTarget" + _), p && h.hasAttribute("morphNormal" + _) === !0 && h.deleteAttribute("morphNormal" + _), r[_] = 0);
      }
      const A = h.morphTargetsRelative ? 1 : 1 - d;
      f.getUniforms().setValue(i, "morphTargetBaseInfluence", A), f.getUniforms().setValue(i, "morphTargetInfluences", r);
    }
  }
  return {
    update: l
  };
}
function Pf(i, e, t, n) {
  let r = /* @__PURE__ */ new WeakMap();
  function s(l) {
    const c = n.render.frame, h = l.geometry, f = e.get(l, h);
    if (r.get(f) !== c && (e.update(f), r.set(f, c)), l.isInstancedMesh && (l.hasEventListener("dispose", a) === !1 && l.addEventListener("dispose", a), r.get(l) !== c && (t.update(l.instanceMatrix, i.ARRAY_BUFFER), l.instanceColor !== null && t.update(l.instanceColor, i.ARRAY_BUFFER), r.set(l, c))), l.isSkinnedMesh) {
      const u = l.skeleton;
      r.get(u) !== c && (u.update(), r.set(u, c));
    }
    return f;
  }
  function o() {
    r = /* @__PURE__ */ new WeakMap();
  }
  function a(l) {
    const c = l.target;
    c.removeEventListener("dispose", a), t.remove(c.instanceMatrix), c.instanceColor !== null && t.remove(c.instanceColor);
  }
  return {
    update: s,
    dispose: o
  };
}
const bo = /* @__PURE__ */ new St(), wo = /* @__PURE__ */ new go(), Ro = /* @__PURE__ */ new gc(), Co = /* @__PURE__ */ new yo(), Ta = [], Aa = [], ba = new Float32Array(16), wa = new Float32Array(9), Ra = new Float32Array(4);
function ri(i, e, t) {
  const n = i[0];
  if (n <= 0 || n > 0)
    return i;
  const r = e * t;
  let s = Ta[r];
  if (s === void 0 && (s = new Float32Array(r), Ta[r] = s), e !== 0) {
    n.toArray(s, 0);
    for (let o = 1, a = 0; o !== e; ++o)
      a += t, i[o].toArray(s, a);
  }
  return s;
}
function rt(i, e) {
  if (i.length !== e.length)
    return !1;
  for (let t = 0, n = i.length; t < n; t++)
    if (i[t] !== e[t])
      return !1;
  return !0;
}
function st(i, e) {
  for (let t = 0, n = e.length; t < n; t++)
    i[t] = e[t];
}
function or(i, e) {
  let t = Aa[e];
  t === void 0 && (t = new Int32Array(e), Aa[e] = t);
  for (let n = 0; n !== e; ++n)
    t[n] = i.allocateTextureUnit();
  return t;
}
function Lf(i, e) {
  const t = this.cache;
  t[0] !== e && (i.uniform1f(this.addr, e), t[0] = e);
}
function Uf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y) && (i.uniform2f(this.addr, e.x, e.y), t[0] = e.x, t[1] = e.y);
  else {
    if (rt(t, e))
      return;
    i.uniform2fv(this.addr, e), st(t, e);
  }
}
function Df(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z) && (i.uniform3f(this.addr, e.x, e.y, e.z), t[0] = e.x, t[1] = e.y, t[2] = e.z);
  else if (e.r !== void 0)
    (t[0] !== e.r || t[1] !== e.g || t[2] !== e.b) && (i.uniform3f(this.addr, e.r, e.g, e.b), t[0] = e.r, t[1] = e.g, t[2] = e.b);
  else {
    if (rt(t, e))
      return;
    i.uniform3fv(this.addr, e), st(t, e);
  }
}
function If(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z || t[3] !== e.w) && (i.uniform4f(this.addr, e.x, e.y, e.z, e.w), t[0] = e.x, t[1] = e.y, t[2] = e.z, t[3] = e.w);
  else {
    if (rt(t, e))
      return;
    i.uniform4fv(this.addr, e), st(t, e);
  }
}
function Nf(i, e) {
  const t = this.cache, n = e.elements;
  if (n === void 0) {
    if (rt(t, e))
      return;
    i.uniformMatrix2fv(this.addr, !1, e), st(t, e);
  } else {
    if (rt(t, n))
      return;
    Ra.set(n), i.uniformMatrix2fv(this.addr, !1, Ra), st(t, n);
  }
}
function Of(i, e) {
  const t = this.cache, n = e.elements;
  if (n === void 0) {
    if (rt(t, e))
      return;
    i.uniformMatrix3fv(this.addr, !1, e), st(t, e);
  } else {
    if (rt(t, n))
      return;
    wa.set(n), i.uniformMatrix3fv(this.addr, !1, wa), st(t, n);
  }
}
function Ff(i, e) {
  const t = this.cache, n = e.elements;
  if (n === void 0) {
    if (rt(t, e))
      return;
    i.uniformMatrix4fv(this.addr, !1, e), st(t, e);
  } else {
    if (rt(t, n))
      return;
    ba.set(n), i.uniformMatrix4fv(this.addr, !1, ba), st(t, n);
  }
}
function Bf(i, e) {
  const t = this.cache;
  t[0] !== e && (i.uniform1i(this.addr, e), t[0] = e);
}
function zf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y) && (i.uniform2i(this.addr, e.x, e.y), t[0] = e.x, t[1] = e.y);
  else {
    if (rt(t, e))
      return;
    i.uniform2iv(this.addr, e), st(t, e);
  }
}
function Hf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z) && (i.uniform3i(this.addr, e.x, e.y, e.z), t[0] = e.x, t[1] = e.y, t[2] = e.z);
  else {
    if (rt(t, e))
      return;
    i.uniform3iv(this.addr, e), st(t, e);
  }
}
function Gf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z || t[3] !== e.w) && (i.uniform4i(this.addr, e.x, e.y, e.z, e.w), t[0] = e.x, t[1] = e.y, t[2] = e.z, t[3] = e.w);
  else {
    if (rt(t, e))
      return;
    i.uniform4iv(this.addr, e), st(t, e);
  }
}
function Vf(i, e) {
  const t = this.cache;
  t[0] !== e && (i.uniform1ui(this.addr, e), t[0] = e);
}
function kf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y) && (i.uniform2ui(this.addr, e.x, e.y), t[0] = e.x, t[1] = e.y);
  else {
    if (rt(t, e))
      return;
    i.uniform2uiv(this.addr, e), st(t, e);
  }
}
function Wf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z) && (i.uniform3ui(this.addr, e.x, e.y, e.z), t[0] = e.x, t[1] = e.y, t[2] = e.z);
  else {
    if (rt(t, e))
      return;
    i.uniform3uiv(this.addr, e), st(t, e);
  }
}
function Xf(i, e) {
  const t = this.cache;
  if (e.x !== void 0)
    (t[0] !== e.x || t[1] !== e.y || t[2] !== e.z || t[3] !== e.w) && (i.uniform4ui(this.addr, e.x, e.y, e.z, e.w), t[0] = e.x, t[1] = e.y, t[2] = e.z, t[3] = e.w);
  else {
    if (rt(t, e))
      return;
    i.uniform4uiv(this.addr, e), st(t, e);
  }
}
function Yf(i, e, t) {
  const n = this.cache, r = t.allocateTextureUnit();
  n[0] !== r && (i.uniform1i(this.addr, r), n[0] = r), t.setTexture2D(e || bo, r);
}
function qf(i, e, t) {
  const n = this.cache, r = t.allocateTextureUnit();
  n[0] !== r && (i.uniform1i(this.addr, r), n[0] = r), t.setTexture3D(e || Ro, r);
}
function jf(i, e, t) {
  const n = this.cache, r = t.allocateTextureUnit();
  n[0] !== r && (i.uniform1i(this.addr, r), n[0] = r), t.setTextureCube(e || Co, r);
}
function Zf(i, e, t) {
  const n = this.cache, r = t.allocateTextureUnit();
  n[0] !== r && (i.uniform1i(this.addr, r), n[0] = r), t.setTexture2DArray(e || wo, r);
}
function Kf(i) {
  switch (i) {
    case 5126:
      return Lf;
    case 35664:
      return Uf;
    case 35665:
      return Df;
    case 35666:
      return If;
    case 35674:
      return Nf;
    case 35675:
      return Of;
    case 35676:
      return Ff;
    case 5124:
    case 35670:
      return Bf;
    case 35667:
    case 35671:
      return zf;
    case 35668:
    case 35672:
      return Hf;
    case 35669:
    case 35673:
      return Gf;
    case 5125:
      return Vf;
    case 36294:
      return kf;
    case 36295:
      return Wf;
    case 36296:
      return Xf;
    case 35678:
    case 36198:
    case 36298:
    case 36306:
    case 35682:
      return Yf;
    case 35679:
    case 36299:
    case 36307:
      return qf;
    case 35680:
    case 36300:
    case 36308:
    case 36293:
      return jf;
    case 36289:
    case 36303:
    case 36311:
    case 36292:
      return Zf;
  }
}
function Jf(i, e) {
  i.uniform1fv(this.addr, e);
}
function $f(i, e) {
  const t = ri(e, this.size, 2);
  i.uniform2fv(this.addr, t);
}
function Qf(i, e) {
  const t = ri(e, this.size, 3);
  i.uniform3fv(this.addr, t);
}
function ed(i, e) {
  const t = ri(e, this.size, 4);
  i.uniform4fv(this.addr, t);
}
function td(i, e) {
  const t = ri(e, this.size, 4);
  i.uniformMatrix2fv(this.addr, !1, t);
}
function nd(i, e) {
  const t = ri(e, this.size, 9);
  i.uniformMatrix3fv(this.addr, !1, t);
}
function id(i, e) {
  const t = ri(e, this.size, 16);
  i.uniformMatrix4fv(this.addr, !1, t);
}
function rd(i, e) {
  i.uniform1iv(this.addr, e);
}
function sd(i, e) {
  i.uniform2iv(this.addr, e);
}
function ad(i, e) {
  i.uniform3iv(this.addr, e);
}
function od(i, e) {
  i.uniform4iv(this.addr, e);
}
function ld(i, e) {
  i.uniform1uiv(this.addr, e);
}
function cd(i, e) {
  i.uniform2uiv(this.addr, e);
}
function hd(i, e) {
  i.uniform3uiv(this.addr, e);
}
function ud(i, e) {
  i.uniform4uiv(this.addr, e);
}
function fd(i, e, t) {
  const n = this.cache, r = e.length, s = or(t, r);
  rt(n, s) || (i.uniform1iv(this.addr, s), st(n, s));
  for (let o = 0; o !== r; ++o)
    t.setTexture2D(e[o] || bo, s[o]);
}
function dd(i, e, t) {
  const n = this.cache, r = e.length, s = or(t, r);
  rt(n, s) || (i.uniform1iv(this.addr, s), st(n, s));
  for (let o = 0; o !== r; ++o)
    t.setTexture3D(e[o] || Ro, s[o]);
}
function pd(i, e, t) {
  const n = this.cache, r = e.length, s = or(t, r);
  rt(n, s) || (i.uniform1iv(this.addr, s), st(n, s));
  for (let o = 0; o !== r; ++o)
    t.setTextureCube(e[o] || Co, s[o]);
}
function md(i, e, t) {
  const n = this.cache, r = e.length, s = or(t, r);
  rt(n, s) || (i.uniform1iv(this.addr, s), st(n, s));
  for (let o = 0; o !== r; ++o)
    t.setTexture2DArray(e[o] || wo, s[o]);
}
function gd(i) {
  switch (i) {
    case 5126:
      return Jf;
    case 35664:
      return $f;
    case 35665:
      return Qf;
    case 35666:
      return ed;
    case 35674:
      return td;
    case 35675:
      return nd;
    case 35676:
      return id;
    case 5124:
    case 35670:
      return rd;
    case 35667:
    case 35671:
      return sd;
    case 35668:
    case 35672:
      return ad;
    case 35669:
    case 35673:
      return od;
    case 5125:
      return ld;
    case 36294:
      return cd;
    case 36295:
      return hd;
    case 36296:
      return ud;
    case 35678:
    case 36198:
    case 36298:
    case 36306:
    case 35682:
      return fd;
    case 35679:
    case 36299:
    case 36307:
      return dd;
    case 35680:
    case 36300:
    case 36308:
    case 36293:
      return pd;
    case 36289:
    case 36303:
    case 36311:
    case 36292:
      return md;
  }
}
class _d {
  constructor(e, t, n) {
    this.id = e, this.addr = n, this.cache = [], this.setValue = Kf(t.type);
  }
}
class vd {
  constructor(e, t, n) {
    this.id = e, this.addr = n, this.cache = [], this.size = t.size, this.setValue = gd(t.type);
  }
}
class xd {
  constructor(e) {
    this.id = e, this.seq = [], this.map = {};
  }
  setValue(e, t, n) {
    const r = this.seq;
    for (let s = 0, o = r.length; s !== o; ++s) {
      const a = r[s];
      a.setValue(e, t[a.id], n);
    }
  }
}
const zr = /(\w+)(\])?(\[|\.)?/g;
function Ca(i, e) {
  i.seq.push(e), i.map[e.id] = e;
}
function Md(i, e, t) {
  const n = i.name, r = n.length;
  for (zr.lastIndex = 0; ; ) {
    const s = zr.exec(n), o = zr.lastIndex;
    let a = s[1];
    const l = s[2] === "]", c = s[3];
    if (l && (a = a | 0), c === void 0 || c === "[" && o + 2 === r) {
      Ca(t, c === void 0 ? new _d(a, i, e) : new vd(a, i, e));
      break;
    } else {
      let f = t.map[a];
      f === void 0 && (f = new xd(a), Ca(t, f)), t = f;
    }
  }
}
class er {
  constructor(e, t) {
    this.seq = [], this.map = {};
    const n = e.getProgramParameter(t, e.ACTIVE_UNIFORMS);
    for (let r = 0; r < n; ++r) {
      const s = e.getActiveUniform(t, r), o = e.getUniformLocation(t, s.name);
      Md(s, o, this);
    }
  }
  setValue(e, t, n, r) {
    const s = this.map[t];
    s !== void 0 && s.setValue(e, n, r);
  }
  setOptional(e, t, n) {
    const r = t[n];
    r !== void 0 && this.setValue(e, n, r);
  }
  static upload(e, t, n, r) {
    for (let s = 0, o = t.length; s !== o; ++s) {
      const a = t[s], l = n[a.id];
      l.needsUpdate !== !1 && a.setValue(e, l.value, r);
    }
  }
  static seqWithValue(e, t) {
    const n = [];
    for (let r = 0, s = e.length; r !== s; ++r) {
      const o = e[r];
      o.id in t && n.push(o);
    }
    return n;
  }
}
function Pa(i, e, t) {
  const n = i.createShader(e);
  return i.shaderSource(n, t), i.compileShader(n), n;
}
let Sd = 0;
function Ed(i, e) {
  const t = i.split(`
`), n = [], r = Math.max(e - 6, 0), s = Math.min(e + 6, t.length);
  for (let o = r; o < s; o++) {
    const a = o + 1;
    n.push(`${a === e ? ">" : " "} ${a}: ${t[o]}`);
  }
  return n.join(`
`);
}
function yd(i) {
  switch (i) {
    case Ft:
      return ["Linear", "( value )"];
    case Oe:
      return ["sRGB", "( value )"];
    default:
      return console.warn("THREE.WebGLProgram: Unsupported color space:", i), ["Linear", "( value )"];
  }
}
function La(i, e, t) {
  const n = i.getShaderParameter(e, i.COMPILE_STATUS), r = i.getShaderInfoLog(e).trim();
  if (n && r === "")
    return "";
  const s = /ERROR: 0:(\d+)/.exec(r);
  if (s) {
    const o = parseInt(s[1]);
    return t.toUpperCase() + `

` + r + `

` + Ed(i.getShaderSource(e), o);
  } else
    return r;
}
function Td(i, e) {
  const t = yd(e);
  return "vec4 " + i + "( vec4 value ) { return LinearTo" + t[0] + t[1] + "; }";
}
function Ad(i, e) {
  let t;
  switch (e) {
    case xl:
      t = "Linear";
      break;
    case Ml:
      t = "Reinhard";
      break;
    case Sl:
      t = "OptimizedCineon";
      break;
    case El:
      t = "ACESFilmic";
      break;
    case yl:
      t = "Custom";
      break;
    default:
      console.warn("THREE.WebGLProgram: Unsupported toneMapping:", e), t = "Linear";
  }
  return "vec3 " + i + "( vec3 color ) { return " + t + "ToneMapping( color ); }";
}
function bd(i) {
  return [
    i.extensionDerivatives || i.envMapCubeUVHeight || i.bumpMap || i.normalMapTangentSpace || i.clearcoatNormalMap || i.flatShading || i.shaderID === "physical" ? "#extension GL_OES_standard_derivatives : enable" : "",
    (i.extensionFragDepth || i.logarithmicDepthBuffer) && i.rendererExtensionFragDepth ? "#extension GL_EXT_frag_depth : enable" : "",
    i.extensionDrawBuffers && i.rendererExtensionDrawBuffers ? "#extension GL_EXT_draw_buffers : require" : "",
    (i.extensionShaderTextureLOD || i.envMap || i.transmission) && i.rendererExtensionShaderTextureLod ? "#extension GL_EXT_shader_texture_lod : enable" : ""
  ].filter(hi).join(`
`);
}
function wd(i) {
  const e = [];
  for (const t in i) {
    const n = i[t];
    n !== !1 && e.push("#define " + t + " " + n);
  }
  return e.join(`
`);
}
function Rd(i, e) {
  const t = {}, n = i.getProgramParameter(e, i.ACTIVE_ATTRIBUTES);
  for (let r = 0; r < n; r++) {
    const s = i.getActiveAttrib(e, r), o = s.name;
    let a = 1;
    s.type === i.FLOAT_MAT2 && (a = 2), s.type === i.FLOAT_MAT3 && (a = 3), s.type === i.FLOAT_MAT4 && (a = 4), t[o] = {
      type: s.type,
      location: i.getAttribLocation(e, o),
      locationSize: a
    };
  }
  return t;
}
function hi(i) {
  return i !== "";
}
function Ua(i, e) {
  const t = e.numSpotLightShadows + e.numSpotLightMaps - e.numSpotLightShadowsWithMaps;
  return i.replace(/NUM_DIR_LIGHTS/g, e.numDirLights).replace(/NUM_SPOT_LIGHTS/g, e.numSpotLights).replace(/NUM_SPOT_LIGHT_MAPS/g, e.numSpotLightMaps).replace(/NUM_SPOT_LIGHT_COORDS/g, t).replace(/NUM_RECT_AREA_LIGHTS/g, e.numRectAreaLights).replace(/NUM_POINT_LIGHTS/g, e.numPointLights).replace(/NUM_HEMI_LIGHTS/g, e.numHemiLights).replace(/NUM_DIR_LIGHT_SHADOWS/g, e.numDirLightShadows).replace(/NUM_SPOT_LIGHT_SHADOWS_WITH_MAPS/g, e.numSpotLightShadowsWithMaps).replace(/NUM_SPOT_LIGHT_SHADOWS/g, e.numSpotLightShadows).replace(/NUM_POINT_LIGHT_SHADOWS/g, e.numPointLightShadows);
}
function Da(i, e) {
  return i.replace(/NUM_CLIPPING_PLANES/g, e.numClippingPlanes).replace(/UNION_CLIPPING_PLANES/g, e.numClippingPlanes - e.numClipIntersection);
}
const Cd = /^[ \t]*#include +<([\w\d./]+)>/gm;
function Qr(i) {
  return i.replace(Cd, Ld);
}
const Pd = /* @__PURE__ */ new Map([
  ["encodings_fragment", "colorspace_fragment"],
  // @deprecated, r154
  ["encodings_pars_fragment", "colorspace_pars_fragment"],
  // @deprecated, r154
  ["output_fragment", "opaque_fragment"]
  // @deprecated, r154
]);
function Ld(i, e) {
  let t = Ne[e];
  if (t === void 0) {
    const n = Pd.get(e);
    if (n !== void 0)
      t = Ne[n], console.warn('THREE.WebGLRenderer: Shader chunk "%s" has been deprecated. Use "%s" instead.', e, n);
    else
      throw new Error("Can not resolve #include <" + e + ">");
  }
  return Qr(t);
}
const Ud = /#pragma unroll_loop_start\s+for\s*\(\s*int\s+i\s*=\s*(\d+)\s*;\s*i\s*<\s*(\d+)\s*;\s*i\s*\+\+\s*\)\s*{([\s\S]+?)}\s+#pragma unroll_loop_end/g;
function Ia(i) {
  return i.replace(Ud, Dd);
}
function Dd(i, e, t, n) {
  let r = "";
  for (let s = parseInt(e); s < parseInt(t); s++)
    r += n.replace(/\[\s*i\s*\]/g, "[ " + s + " ]").replace(/UNROLLED_LOOP_INDEX/g, s);
  return r;
}
function Na(i) {
  let e = "precision " + i.precision + ` float;
precision ` + i.precision + " int;";
  return i.precision === "highp" ? e += `
#define HIGH_PRECISION` : i.precision === "mediump" ? e += `
#define MEDIUM_PRECISION` : i.precision === "lowp" && (e += `
#define LOW_PRECISION`), e;
}
function Id(i) {
  let e = "SHADOWMAP_TYPE_BASIC";
  return i.shadowMapType === $a ? e = "SHADOWMAP_TYPE_PCF" : i.shadowMapType === Jo ? e = "SHADOWMAP_TYPE_PCF_SOFT" : i.shadowMapType === Yt && (e = "SHADOWMAP_TYPE_VSM"), e;
}
function Nd(i) {
  let e = "ENVMAP_TYPE_CUBE";
  if (i.envMap)
    switch (i.envMapMode) {
      case ei:
      case ti:
        e = "ENVMAP_TYPE_CUBE";
        break;
      case rr:
        e = "ENVMAP_TYPE_CUBE_UV";
        break;
    }
  return e;
}
function Od(i) {
  let e = "ENVMAP_MODE_REFLECTION";
  if (i.envMap)
    switch (i.envMapMode) {
      case ti:
        e = "ENVMAP_MODE_REFRACTION";
        break;
    }
  return e;
}
function Fd(i) {
  let e = "ENVMAP_BLENDING_NONE";
  if (i.envMap)
    switch (i.combine) {
      case to:
        e = "ENVMAP_BLENDING_MULTIPLY";
        break;
      case _l:
        e = "ENVMAP_BLENDING_MIX";
        break;
      case vl:
        e = "ENVMAP_BLENDING_ADD";
        break;
    }
  return e;
}
function Bd(i) {
  const e = i.envMapCubeUVHeight;
  if (e === null)
    return null;
  const t = Math.log2(e) - 2, n = 1 / e;
  return { texelWidth: 1 / (3 * Math.max(Math.pow(2, t), 7 * 16)), texelHeight: n, maxMip: t };
}
function zd(i, e, t, n) {
  const r = i.getContext(), s = t.defines;
  let o = t.vertexShader, a = t.fragmentShader;
  const l = Id(t), c = Nd(t), h = Od(t), f = Fd(t), u = Bd(t), m = t.isWebGL2 ? "" : bd(t), g = wd(s), x = r.createProgram();
  let p, d, A = t.glslVersion ? "#version " + t.glslVersion + `
` : "";
  t.isRawShaderMaterial ? (p = [
    "#define SHADER_TYPE " + t.shaderType,
    "#define SHADER_NAME " + t.shaderName,
    g
  ].filter(hi).join(`
`), p.length > 0 && (p += `
`), d = [
    m,
    "#define SHADER_TYPE " + t.shaderType,
    "#define SHADER_NAME " + t.shaderName,
    g
  ].filter(hi).join(`
`), d.length > 0 && (d += `
`)) : (p = [
    Na(t),
    "#define SHADER_TYPE " + t.shaderType,
    "#define SHADER_NAME " + t.shaderName,
    g,
    t.instancing ? "#define USE_INSTANCING" : "",
    t.instancingColor ? "#define USE_INSTANCING_COLOR" : "",
    t.useFog && t.fog ? "#define USE_FOG" : "",
    t.useFog && t.fogExp2 ? "#define FOG_EXP2" : "",
    t.map ? "#define USE_MAP" : "",
    t.envMap ? "#define USE_ENVMAP" : "",
    t.envMap ? "#define " + h : "",
    t.lightMap ? "#define USE_LIGHTMAP" : "",
    t.aoMap ? "#define USE_AOMAP" : "",
    t.bumpMap ? "#define USE_BUMPMAP" : "",
    t.normalMap ? "#define USE_NORMALMAP" : "",
    t.normalMapObjectSpace ? "#define USE_NORMALMAP_OBJECTSPACE" : "",
    t.normalMapTangentSpace ? "#define USE_NORMALMAP_TANGENTSPACE" : "",
    t.displacementMap ? "#define USE_DISPLACEMENTMAP" : "",
    t.emissiveMap ? "#define USE_EMISSIVEMAP" : "",
    t.anisotropyMap ? "#define USE_ANISOTROPYMAP" : "",
    t.clearcoatMap ? "#define USE_CLEARCOATMAP" : "",
    t.clearcoatRoughnessMap ? "#define USE_CLEARCOAT_ROUGHNESSMAP" : "",
    t.clearcoatNormalMap ? "#define USE_CLEARCOAT_NORMALMAP" : "",
    t.iridescenceMap ? "#define USE_IRIDESCENCEMAP" : "",
    t.iridescenceThicknessMap ? "#define USE_IRIDESCENCE_THICKNESSMAP" : "",
    t.specularMap ? "#define USE_SPECULARMAP" : "",
    t.specularColorMap ? "#define USE_SPECULAR_COLORMAP" : "",
    t.specularIntensityMap ? "#define USE_SPECULAR_INTENSITYMAP" : "",
    t.roughnessMap ? "#define USE_ROUGHNESSMAP" : "",
    t.metalnessMap ? "#define USE_METALNESSMAP" : "",
    t.alphaMap ? "#define USE_ALPHAMAP" : "",
    t.alphaHash ? "#define USE_ALPHAHASH" : "",
    t.transmission ? "#define USE_TRANSMISSION" : "",
    t.transmissionMap ? "#define USE_TRANSMISSIONMAP" : "",
    t.thicknessMap ? "#define USE_THICKNESSMAP" : "",
    t.sheenColorMap ? "#define USE_SHEEN_COLORMAP" : "",
    t.sheenRoughnessMap ? "#define USE_SHEEN_ROUGHNESSMAP" : "",
    //
    t.mapUv ? "#define MAP_UV " + t.mapUv : "",
    t.alphaMapUv ? "#define ALPHAMAP_UV " + t.alphaMapUv : "",
    t.lightMapUv ? "#define LIGHTMAP_UV " + t.lightMapUv : "",
    t.aoMapUv ? "#define AOMAP_UV " + t.aoMapUv : "",
    t.emissiveMapUv ? "#define EMISSIVEMAP_UV " + t.emissiveMapUv : "",
    t.bumpMapUv ? "#define BUMPMAP_UV " + t.bumpMapUv : "",
    t.normalMapUv ? "#define NORMALMAP_UV " + t.normalMapUv : "",
    t.displacementMapUv ? "#define DISPLACEMENTMAP_UV " + t.displacementMapUv : "",
    t.metalnessMapUv ? "#define METALNESSMAP_UV " + t.metalnessMapUv : "",
    t.roughnessMapUv ? "#define ROUGHNESSMAP_UV " + t.roughnessMapUv : "",
    t.anisotropyMapUv ? "#define ANISOTROPYMAP_UV " + t.anisotropyMapUv : "",
    t.clearcoatMapUv ? "#define CLEARCOATMAP_UV " + t.clearcoatMapUv : "",
    t.clearcoatNormalMapUv ? "#define CLEARCOAT_NORMALMAP_UV " + t.clearcoatNormalMapUv : "",
    t.clearcoatRoughnessMapUv ? "#define CLEARCOAT_ROUGHNESSMAP_UV " + t.clearcoatRoughnessMapUv : "",
    t.iridescenceMapUv ? "#define IRIDESCENCEMAP_UV " + t.iridescenceMapUv : "",
    t.iridescenceThicknessMapUv ? "#define IRIDESCENCE_THICKNESSMAP_UV " + t.iridescenceThicknessMapUv : "",
    t.sheenColorMapUv ? "#define SHEEN_COLORMAP_UV " + t.sheenColorMapUv : "",
    t.sheenRoughnessMapUv ? "#define SHEEN_ROUGHNESSMAP_UV " + t.sheenRoughnessMapUv : "",
    t.specularMapUv ? "#define SPECULARMAP_UV " + t.specularMapUv : "",
    t.specularColorMapUv ? "#define SPECULAR_COLORMAP_UV " + t.specularColorMapUv : "",
    t.specularIntensityMapUv ? "#define SPECULAR_INTENSITYMAP_UV " + t.specularIntensityMapUv : "",
    t.transmissionMapUv ? "#define TRANSMISSIONMAP_UV " + t.transmissionMapUv : "",
    t.thicknessMapUv ? "#define THICKNESSMAP_UV " + t.thicknessMapUv : "",
    //
    t.vertexTangents && t.flatShading === !1 ? "#define USE_TANGENT" : "",
    t.vertexColors ? "#define USE_COLOR" : "",
    t.vertexAlphas ? "#define USE_COLOR_ALPHA" : "",
    t.vertexUv1s ? "#define USE_UV1" : "",
    t.vertexUv2s ? "#define USE_UV2" : "",
    t.vertexUv3s ? "#define USE_UV3" : "",
    t.pointsUvs ? "#define USE_POINTS_UV" : "",
    t.flatShading ? "#define FLAT_SHADED" : "",
    t.skinning ? "#define USE_SKINNING" : "",
    t.morphTargets ? "#define USE_MORPHTARGETS" : "",
    t.morphNormals && t.flatShading === !1 ? "#define USE_MORPHNORMALS" : "",
    t.morphColors && t.isWebGL2 ? "#define USE_MORPHCOLORS" : "",
    t.morphTargetsCount > 0 && t.isWebGL2 ? "#define MORPHTARGETS_TEXTURE" : "",
    t.morphTargetsCount > 0 && t.isWebGL2 ? "#define MORPHTARGETS_TEXTURE_STRIDE " + t.morphTextureStride : "",
    t.morphTargetsCount > 0 && t.isWebGL2 ? "#define MORPHTARGETS_COUNT " + t.morphTargetsCount : "",
    t.doubleSided ? "#define DOUBLE_SIDED" : "",
    t.flipSided ? "#define FLIP_SIDED" : "",
    t.shadowMapEnabled ? "#define USE_SHADOWMAP" : "",
    t.shadowMapEnabled ? "#define " + l : "",
    t.sizeAttenuation ? "#define USE_SIZEATTENUATION" : "",
    t.useLegacyLights ? "#define LEGACY_LIGHTS" : "",
    t.logarithmicDepthBuffer ? "#define USE_LOGDEPTHBUF" : "",
    t.logarithmicDepthBuffer && t.rendererExtensionFragDepth ? "#define USE_LOGDEPTHBUF_EXT" : "",
    "uniform mat4 modelMatrix;",
    "uniform mat4 modelViewMatrix;",
    "uniform mat4 projectionMatrix;",
    "uniform mat4 viewMatrix;",
    "uniform mat3 normalMatrix;",
    "uniform vec3 cameraPosition;",
    "uniform bool isOrthographic;",
    "#ifdef USE_INSTANCING",
    "	attribute mat4 instanceMatrix;",
    "#endif",
    "#ifdef USE_INSTANCING_COLOR",
    "	attribute vec3 instanceColor;",
    "#endif",
    "attribute vec3 position;",
    "attribute vec3 normal;",
    "attribute vec2 uv;",
    "#ifdef USE_UV1",
    "	attribute vec2 uv1;",
    "#endif",
    "#ifdef USE_UV2",
    "	attribute vec2 uv2;",
    "#endif",
    "#ifdef USE_UV3",
    "	attribute vec2 uv3;",
    "#endif",
    "#ifdef USE_TANGENT",
    "	attribute vec4 tangent;",
    "#endif",
    "#if defined( USE_COLOR_ALPHA )",
    "	attribute vec4 color;",
    "#elif defined( USE_COLOR )",
    "	attribute vec3 color;",
    "#endif",
    "#if ( defined( USE_MORPHTARGETS ) && ! defined( MORPHTARGETS_TEXTURE ) )",
    "	attribute vec3 morphTarget0;",
    "	attribute vec3 morphTarget1;",
    "	attribute vec3 morphTarget2;",
    "	attribute vec3 morphTarget3;",
    "	#ifdef USE_MORPHNORMALS",
    "		attribute vec3 morphNormal0;",
    "		attribute vec3 morphNormal1;",
    "		attribute vec3 morphNormal2;",
    "		attribute vec3 morphNormal3;",
    "	#else",
    "		attribute vec3 morphTarget4;",
    "		attribute vec3 morphTarget5;",
    "		attribute vec3 morphTarget6;",
    "		attribute vec3 morphTarget7;",
    "	#endif",
    "#endif",
    "#ifdef USE_SKINNING",
    "	attribute vec4 skinIndex;",
    "	attribute vec4 skinWeight;",
    "#endif",
    `
`
  ].filter(hi).join(`
`), d = [
    m,
    Na(t),
    "#define SHADER_TYPE " + t.shaderType,
    "#define SHADER_NAME " + t.shaderName,
    g,
    t.useFog && t.fog ? "#define USE_FOG" : "",
    t.useFog && t.fogExp2 ? "#define FOG_EXP2" : "",
    t.map ? "#define USE_MAP" : "",
    t.matcap ? "#define USE_MATCAP" : "",
    t.envMap ? "#define USE_ENVMAP" : "",
    t.envMap ? "#define " + c : "",
    t.envMap ? "#define " + h : "",
    t.envMap ? "#define " + f : "",
    u ? "#define CUBEUV_TEXEL_WIDTH " + u.texelWidth : "",
    u ? "#define CUBEUV_TEXEL_HEIGHT " + u.texelHeight : "",
    u ? "#define CUBEUV_MAX_MIP " + u.maxMip + ".0" : "",
    t.lightMap ? "#define USE_LIGHTMAP" : "",
    t.aoMap ? "#define USE_AOMAP" : "",
    t.bumpMap ? "#define USE_BUMPMAP" : "",
    t.normalMap ? "#define USE_NORMALMAP" : "",
    t.normalMapObjectSpace ? "#define USE_NORMALMAP_OBJECTSPACE" : "",
    t.normalMapTangentSpace ? "#define USE_NORMALMAP_TANGENTSPACE" : "",
    t.emissiveMap ? "#define USE_EMISSIVEMAP" : "",
    t.anisotropy ? "#define USE_ANISOTROPY" : "",
    t.anisotropyMap ? "#define USE_ANISOTROPYMAP" : "",
    t.clearcoat ? "#define USE_CLEARCOAT" : "",
    t.clearcoatMap ? "#define USE_CLEARCOATMAP" : "",
    t.clearcoatRoughnessMap ? "#define USE_CLEARCOAT_ROUGHNESSMAP" : "",
    t.clearcoatNormalMap ? "#define USE_CLEARCOAT_NORMALMAP" : "",
    t.iridescence ? "#define USE_IRIDESCENCE" : "",
    t.iridescenceMap ? "#define USE_IRIDESCENCEMAP" : "",
    t.iridescenceThicknessMap ? "#define USE_IRIDESCENCE_THICKNESSMAP" : "",
    t.specularMap ? "#define USE_SPECULARMAP" : "",
    t.specularColorMap ? "#define USE_SPECULAR_COLORMAP" : "",
    t.specularIntensityMap ? "#define USE_SPECULAR_INTENSITYMAP" : "",
    t.roughnessMap ? "#define USE_ROUGHNESSMAP" : "",
    t.metalnessMap ? "#define USE_METALNESSMAP" : "",
    t.alphaMap ? "#define USE_ALPHAMAP" : "",
    t.alphaTest ? "#define USE_ALPHATEST" : "",
    t.alphaHash ? "#define USE_ALPHAHASH" : "",
    t.sheen ? "#define USE_SHEEN" : "",
    t.sheenColorMap ? "#define USE_SHEEN_COLORMAP" : "",
    t.sheenRoughnessMap ? "#define USE_SHEEN_ROUGHNESSMAP" : "",
    t.transmission ? "#define USE_TRANSMISSION" : "",
    t.transmissionMap ? "#define USE_TRANSMISSIONMAP" : "",
    t.thicknessMap ? "#define USE_THICKNESSMAP" : "",
    t.vertexTangents && t.flatShading === !1 ? "#define USE_TANGENT" : "",
    t.vertexColors || t.instancingColor ? "#define USE_COLOR" : "",
    t.vertexAlphas ? "#define USE_COLOR_ALPHA" : "",
    t.vertexUv1s ? "#define USE_UV1" : "",
    t.vertexUv2s ? "#define USE_UV2" : "",
    t.vertexUv3s ? "#define USE_UV3" : "",
    t.pointsUvs ? "#define USE_POINTS_UV" : "",
    t.gradientMap ? "#define USE_GRADIENTMAP" : "",
    t.flatShading ? "#define FLAT_SHADED" : "",
    t.doubleSided ? "#define DOUBLE_SIDED" : "",
    t.flipSided ? "#define FLIP_SIDED" : "",
    t.shadowMapEnabled ? "#define USE_SHADOWMAP" : "",
    t.shadowMapEnabled ? "#define " + l : "",
    t.premultipliedAlpha ? "#define PREMULTIPLIED_ALPHA" : "",
    t.useLegacyLights ? "#define LEGACY_LIGHTS" : "",
    t.logarithmicDepthBuffer ? "#define USE_LOGDEPTHBUF" : "",
    t.logarithmicDepthBuffer && t.rendererExtensionFragDepth ? "#define USE_LOGDEPTHBUF_EXT" : "",
    "uniform mat4 viewMatrix;",
    "uniform vec3 cameraPosition;",
    "uniform bool isOrthographic;",
    t.toneMapping !== an ? "#define TONE_MAPPING" : "",
    t.toneMapping !== an ? Ne.tonemapping_pars_fragment : "",
    // this code is required here because it is used by the toneMapping() function defined below
    t.toneMapping !== an ? Ad("toneMapping", t.toneMapping) : "",
    t.dithering ? "#define DITHERING" : "",
    t.opaque ? "#define OPAQUE" : "",
    Ne.colorspace_pars_fragment,
    // this code is required here because it is used by the various encoding/decoding function defined below
    Td("linearToOutputTexel", t.outputColorSpace),
    t.useDepthPacking ? "#define DEPTH_PACKING " + t.depthPacking : "",
    `
`
  ].filter(hi).join(`
`)), o = Qr(o), o = Ua(o, t), o = Da(o, t), a = Qr(a), a = Ua(a, t), a = Da(a, t), o = Ia(o), a = Ia(a), t.isWebGL2 && t.isRawShaderMaterial !== !0 && (A = `#version 300 es
`, p = [
    "precision mediump sampler2DArray;",
    "#define attribute in",
    "#define varying out",
    "#define texture2D texture"
  ].join(`
`) + `
` + p, d = [
    "#define varying in",
    t.glslVersion === ea ? "" : "layout(location = 0) out highp vec4 pc_fragColor;",
    t.glslVersion === ea ? "" : "#define gl_FragColor pc_fragColor",
    "#define gl_FragDepthEXT gl_FragDepth",
    "#define texture2D texture",
    "#define textureCube texture",
    "#define texture2DProj textureProj",
    "#define texture2DLodEXT textureLod",
    "#define texture2DProjLodEXT textureProjLod",
    "#define textureCubeLodEXT textureLod",
    "#define texture2DGradEXT textureGrad",
    "#define texture2DProjGradEXT textureProjGrad",
    "#define textureCubeGradEXT textureGrad"
  ].join(`
`) + `
` + d);
  const _ = A + p + o, T = A + d + a, C = Pa(r, r.VERTEX_SHADER, _), L = Pa(r, r.FRAGMENT_SHADER, T);
  if (r.attachShader(x, C), r.attachShader(x, L), t.index0AttributeName !== void 0 ? r.bindAttribLocation(x, 0, t.index0AttributeName) : t.morphTargets === !0 && r.bindAttribLocation(x, 0, "position"), r.linkProgram(x), i.debug.checkShaderErrors) {
    const M = r.getProgramInfoLog(x).trim(), y = r.getShaderInfoLog(C).trim(), Y = r.getShaderInfoLog(L).trim();
    let ce = !0, B = !0;
    if (r.getProgramParameter(x, r.LINK_STATUS) === !1)
      if (ce = !1, typeof i.debug.onShaderError == "function")
        i.debug.onShaderError(r, x, C, L);
      else {
        const H = La(r, C, "vertex"), G = La(r, L, "fragment");
        console.error(
          "THREE.WebGLProgram: Shader Error " + r.getError() + " - VALIDATE_STATUS " + r.getProgramParameter(x, r.VALIDATE_STATUS) + `

Program Info Log: ` + M + `
` + H + `
` + G
        );
      }
    else
      M !== "" ? console.warn("THREE.WebGLProgram: Program Info Log:", M) : (y === "" || Y === "") && (B = !1);
    B && (this.diagnostics = {
      runnable: ce,
      programLog: M,
      vertexShader: {
        log: y,
        prefix: p
      },
      fragmentShader: {
        log: Y,
        prefix: d
      }
    });
  }
  r.deleteShader(C), r.deleteShader(L);
  let w;
  this.getUniforms = function() {
    return w === void 0 && (w = new er(r, x)), w;
  };
  let V;
  return this.getAttributes = function() {
    return V === void 0 && (V = Rd(r, x)), V;
  }, this.destroy = function() {
    n.releaseStatesOfProgram(this), r.deleteProgram(x), this.program = void 0;
  }, this.type = t.shaderType, this.name = t.shaderName, this.id = Sd++, this.cacheKey = e, this.usedTimes = 1, this.program = x, this.vertexShader = C, this.fragmentShader = L, this;
}
let Hd = 0;
class Gd {
  constructor() {
    this.shaderCache = /* @__PURE__ */ new Map(), this.materialCache = /* @__PURE__ */ new Map();
  }
  update(e) {
    const t = e.vertexShader, n = e.fragmentShader, r = this._getShaderStage(t), s = this._getShaderStage(n), o = this._getShaderCacheForMaterial(e);
    return o.has(r) === !1 && (o.add(r), r.usedTimes++), o.has(s) === !1 && (o.add(s), s.usedTimes++), this;
  }
  remove(e) {
    const t = this.materialCache.get(e);
    for (const n of t)
      n.usedTimes--, n.usedTimes === 0 && this.shaderCache.delete(n.code);
    return this.materialCache.delete(e), this;
  }
  getVertexShaderID(e) {
    return this._getShaderStage(e.vertexShader).id;
  }
  getFragmentShaderID(e) {
    return this._getShaderStage(e.fragmentShader).id;
  }
  dispose() {
    this.shaderCache.clear(), this.materialCache.clear();
  }
  _getShaderCacheForMaterial(e) {
    const t = this.materialCache;
    let n = t.get(e);
    return n === void 0 && (n = /* @__PURE__ */ new Set(), t.set(e, n)), n;
  }
  _getShaderStage(e) {
    const t = this.shaderCache;
    let n = t.get(e);
    return n === void 0 && (n = new Vd(e), t.set(e, n)), n;
  }
}
class Vd {
  constructor(e) {
    this.id = Hd++, this.code = e, this.usedTimes = 0;
  }
}
function kd(i, e, t, n, r, s, o) {
  const a = new hs(), l = new Gd(), c = [], h = r.isWebGL2, f = r.logarithmicDepthBuffer, u = r.vertexTextures;
  let m = r.precision;
  const g = {
    MeshDepthMaterial: "depth",
    MeshDistanceMaterial: "distanceRGBA",
    MeshNormalMaterial: "normal",
    MeshBasicMaterial: "basic",
    MeshLambertMaterial: "lambert",
    MeshPhongMaterial: "phong",
    MeshToonMaterial: "toon",
    MeshStandardMaterial: "physical",
    MeshPhysicalMaterial: "physical",
    MeshMatcapMaterial: "matcap",
    LineBasicMaterial: "basic",
    LineDashedMaterial: "dashed",
    PointsMaterial: "points",
    ShadowMaterial: "shadow",
    SpriteMaterial: "sprite"
  };
  function x(M) {
    return M === 0 ? "uv" : `uv${M}`;
  }
  function p(M, y, Y, ce, B) {
    const H = ce.fog, G = B.geometry, Q = M.isMeshStandardMaterial ? ce.environment : null, X = (M.isMeshStandardMaterial ? t : e).get(M.envMap || Q), j = X && X.mapping === rr ? X.image.height : null, J = g[M.type];
    M.precision !== null && (m = r.getMaxPrecision(M.precision), m !== M.precision && console.warn("THREE.WebGLProgram.getParameters:", M.precision, "not supported, using", m, "instead."));
    const ee = G.morphAttributes.position || G.morphAttributes.normal || G.morphAttributes.color, I = ee !== void 0 ? ee.length : 0;
    let q = 0;
    G.morphAttributes.position !== void 0 && (q = 1), G.morphAttributes.normal !== void 0 && (q = 2), G.morphAttributes.color !== void 0 && (q = 3);
    let pe, me, ve, Ae;
    if (J) {
      const je = It[J];
      pe = je.vertexShader, me = je.fragmentShader;
    } else
      pe = M.vertexShader, me = M.fragmentShader, l.update(M), ve = l.getVertexShaderID(M), Ae = l.getFragmentShaderID(M);
    const be = i.getRenderTarget(), Te = B.isInstancedMesh === !0, ke = !!M.map, Ye = !!M.matcap, we = !!X, b = !!M.aoMap, le = !!M.lightMap, Z = !!M.bumpMap, re = !!M.normalMap, $ = !!M.displacementMap, ye = !!M.emissiveMap, xe = !!M.metalnessMap, Me = !!M.roughnessMap, Ce = M.anisotropy > 0, ze = M.clearcoat > 0, Ze = M.iridescence > 0, E = M.sheen > 0, v = M.transmission > 0, O = Ce && !!M.anisotropyMap, ie = ze && !!M.clearcoatMap, te = ze && !!M.clearcoatNormalMap, se = ze && !!M.clearcoatRoughnessMap, Ee = Ze && !!M.iridescenceMap, ae = Ze && !!M.iridescenceThicknessMap, z = E && !!M.sheenColorMap, R = E && !!M.sheenRoughnessMap, K = !!M.specularMap, _e = !!M.specularColorMap, he = !!M.specularIntensityMap, ge = v && !!M.transmissionMap, De = v && !!M.thicknessMap, Ve = !!M.gradientMap, P = !!M.alphaMap, fe = M.alphaTest > 0, F = !!M.alphaHash, ne = !!M.extensions, de = !!G.attributes.uv1, Fe = !!G.attributes.uv2, Xe = !!G.attributes.uv3;
    let qe = an;
    return M.toneMapped && (be === null || be.isXRRenderTarget === !0) && (qe = i.toneMapping), {
      isWebGL2: h,
      shaderID: J,
      shaderType: M.type,
      shaderName: M.name,
      vertexShader: pe,
      fragmentShader: me,
      defines: M.defines,
      customVertexShaderID: ve,
      customFragmentShaderID: Ae,
      isRawShaderMaterial: M.isRawShaderMaterial === !0,
      glslVersion: M.glslVersion,
      precision: m,
      instancing: Te,
      instancingColor: Te && B.instanceColor !== null,
      supportsVertexTextures: u,
      outputColorSpace: be === null ? i.outputColorSpace : be.isXRRenderTarget === !0 ? be.texture.colorSpace : Ft,
      map: ke,
      matcap: Ye,
      envMap: we,
      envMapMode: we && X.mapping,
      envMapCubeUVHeight: j,
      aoMap: b,
      lightMap: le,
      bumpMap: Z,
      normalMap: re,
      displacementMap: u && $,
      emissiveMap: ye,
      normalMapObjectSpace: re && M.normalMapType === Ol,
      normalMapTangentSpace: re && M.normalMapType === ho,
      metalnessMap: xe,
      roughnessMap: Me,
      anisotropy: Ce,
      anisotropyMap: O,
      clearcoat: ze,
      clearcoatMap: ie,
      clearcoatNormalMap: te,
      clearcoatRoughnessMap: se,
      iridescence: Ze,
      iridescenceMap: Ee,
      iridescenceThicknessMap: ae,
      sheen: E,
      sheenColorMap: z,
      sheenRoughnessMap: R,
      specularMap: K,
      specularColorMap: _e,
      specularIntensityMap: he,
      transmission: v,
      transmissionMap: ge,
      thicknessMap: De,
      gradientMap: Ve,
      opaque: M.transparent === !1 && M.blending === $n,
      alphaMap: P,
      alphaTest: fe,
      alphaHash: F,
      combine: M.combine,
      //
      mapUv: ke && x(M.map.channel),
      aoMapUv: b && x(M.aoMap.channel),
      lightMapUv: le && x(M.lightMap.channel),
      bumpMapUv: Z && x(M.bumpMap.channel),
      normalMapUv: re && x(M.normalMap.channel),
      displacementMapUv: $ && x(M.displacementMap.channel),
      emissiveMapUv: ye && x(M.emissiveMap.channel),
      metalnessMapUv: xe && x(M.metalnessMap.channel),
      roughnessMapUv: Me && x(M.roughnessMap.channel),
      anisotropyMapUv: O && x(M.anisotropyMap.channel),
      clearcoatMapUv: ie && x(M.clearcoatMap.channel),
      clearcoatNormalMapUv: te && x(M.clearcoatNormalMap.channel),
      clearcoatRoughnessMapUv: se && x(M.clearcoatRoughnessMap.channel),
      iridescenceMapUv: Ee && x(M.iridescenceMap.channel),
      iridescenceThicknessMapUv: ae && x(M.iridescenceThicknessMap.channel),
      sheenColorMapUv: z && x(M.sheenColorMap.channel),
      sheenRoughnessMapUv: R && x(M.sheenRoughnessMap.channel),
      specularMapUv: K && x(M.specularMap.channel),
      specularColorMapUv: _e && x(M.specularColorMap.channel),
      specularIntensityMapUv: he && x(M.specularIntensityMap.channel),
      transmissionMapUv: ge && x(M.transmissionMap.channel),
      thicknessMapUv: De && x(M.thicknessMap.channel),
      alphaMapUv: P && x(M.alphaMap.channel),
      //
      vertexTangents: !!G.attributes.tangent && (re || Ce),
      vertexColors: M.vertexColors,
      vertexAlphas: M.vertexColors === !0 && !!G.attributes.color && G.attributes.color.itemSize === 4,
      vertexUv1s: de,
      vertexUv2s: Fe,
      vertexUv3s: Xe,
      pointsUvs: B.isPoints === !0 && !!G.attributes.uv && (ke || P),
      fog: !!H,
      useFog: M.fog === !0,
      fogExp2: H && H.isFogExp2,
      flatShading: M.flatShading === !0,
      sizeAttenuation: M.sizeAttenuation === !0,
      logarithmicDepthBuffer: f,
      skinning: B.isSkinnedMesh === !0,
      morphTargets: G.morphAttributes.position !== void 0,
      morphNormals: G.morphAttributes.normal !== void 0,
      morphColors: G.morphAttributes.color !== void 0,
      morphTargetsCount: I,
      morphTextureStride: q,
      numDirLights: y.directional.length,
      numPointLights: y.point.length,
      numSpotLights: y.spot.length,
      numSpotLightMaps: y.spotLightMap.length,
      numRectAreaLights: y.rectArea.length,
      numHemiLights: y.hemi.length,
      numDirLightShadows: y.directionalShadowMap.length,
      numPointLightShadows: y.pointShadowMap.length,
      numSpotLightShadows: y.spotShadowMap.length,
      numSpotLightShadowsWithMaps: y.numSpotLightShadowsWithMaps,
      numClippingPlanes: o.numPlanes,
      numClipIntersection: o.numIntersection,
      dithering: M.dithering,
      shadowMapEnabled: i.shadowMap.enabled && Y.length > 0,
      shadowMapType: i.shadowMap.type,
      toneMapping: qe,
      useLegacyLights: i._useLegacyLights,
      premultipliedAlpha: M.premultipliedAlpha,
      doubleSided: M.side === Mt,
      flipSided: M.side === gt,
      useDepthPacking: M.depthPacking >= 0,
      depthPacking: M.depthPacking || 0,
      index0AttributeName: M.index0AttributeName,
      extensionDerivatives: ne && M.extensions.derivatives === !0,
      extensionFragDepth: ne && M.extensions.fragDepth === !0,
      extensionDrawBuffers: ne && M.extensions.drawBuffers === !0,
      extensionShaderTextureLOD: ne && M.extensions.shaderTextureLOD === !0,
      rendererExtensionFragDepth: h || n.has("EXT_frag_depth"),
      rendererExtensionDrawBuffers: h || n.has("WEBGL_draw_buffers"),
      rendererExtensionShaderTextureLod: h || n.has("EXT_shader_texture_lod"),
      customProgramCacheKey: M.customProgramCacheKey()
    };
  }
  function d(M) {
    const y = [];
    if (M.shaderID ? y.push(M.shaderID) : (y.push(M.customVertexShaderID), y.push(M.customFragmentShaderID)), M.defines !== void 0)
      for (const Y in M.defines)
        y.push(Y), y.push(M.defines[Y]);
    return M.isRawShaderMaterial === !1 && (A(y, M), _(y, M), y.push(i.outputColorSpace)), y.push(M.customProgramCacheKey), y.join();
  }
  function A(M, y) {
    M.push(y.precision), M.push(y.outputColorSpace), M.push(y.envMapMode), M.push(y.envMapCubeUVHeight), M.push(y.mapUv), M.push(y.alphaMapUv), M.push(y.lightMapUv), M.push(y.aoMapUv), M.push(y.bumpMapUv), M.push(y.normalMapUv), M.push(y.displacementMapUv), M.push(y.emissiveMapUv), M.push(y.metalnessMapUv), M.push(y.roughnessMapUv), M.push(y.anisotropyMapUv), M.push(y.clearcoatMapUv), M.push(y.clearcoatNormalMapUv), M.push(y.clearcoatRoughnessMapUv), M.push(y.iridescenceMapUv), M.push(y.iridescenceThicknessMapUv), M.push(y.sheenColorMapUv), M.push(y.sheenRoughnessMapUv), M.push(y.specularMapUv), M.push(y.specularColorMapUv), M.push(y.specularIntensityMapUv), M.push(y.transmissionMapUv), M.push(y.thicknessMapUv), M.push(y.combine), M.push(y.fogExp2), M.push(y.sizeAttenuation), M.push(y.morphTargetsCount), M.push(y.morphAttributeCount), M.push(y.numDirLights), M.push(y.numPointLights), M.push(y.numSpotLights), M.push(y.numSpotLightMaps), M.push(y.numHemiLights), M.push(y.numRectAreaLights), M.push(y.numDirLightShadows), M.push(y.numPointLightShadows), M.push(y.numSpotLightShadows), M.push(y.numSpotLightShadowsWithMaps), M.push(y.shadowMapType), M.push(y.toneMapping), M.push(y.numClippingPlanes), M.push(y.numClipIntersection), M.push(y.depthPacking);
  }
  function _(M, y) {
    a.disableAll(), y.isWebGL2 && a.enable(0), y.supportsVertexTextures && a.enable(1), y.instancing && a.enable(2), y.instancingColor && a.enable(3), y.matcap && a.enable(4), y.envMap && a.enable(5), y.normalMapObjectSpace && a.enable(6), y.normalMapTangentSpace && a.enable(7), y.clearcoat && a.enable(8), y.iridescence && a.enable(9), y.alphaTest && a.enable(10), y.vertexColors && a.enable(11), y.vertexAlphas && a.enable(12), y.vertexUv1s && a.enable(13), y.vertexUv2s && a.enable(14), y.vertexUv3s && a.enable(15), y.vertexTangents && a.enable(16), y.anisotropy && a.enable(17), M.push(a.mask), a.disableAll(), y.fog && a.enable(0), y.useFog && a.enable(1), y.flatShading && a.enable(2), y.logarithmicDepthBuffer && a.enable(3), y.skinning && a.enable(4), y.morphTargets && a.enable(5), y.morphNormals && a.enable(6), y.morphColors && a.enable(7), y.premultipliedAlpha && a.enable(8), y.shadowMapEnabled && a.enable(9), y.useLegacyLights && a.enable(10), y.doubleSided && a.enable(11), y.flipSided && a.enable(12), y.useDepthPacking && a.enable(13), y.dithering && a.enable(14), y.transmission && a.enable(15), y.sheen && a.enable(16), y.opaque && a.enable(17), y.pointsUvs && a.enable(18), M.push(a.mask);
  }
  function T(M) {
    const y = g[M.type];
    let Y;
    if (y) {
      const ce = It[y];
      Y = Rc.clone(ce.uniforms);
    } else
      Y = M.uniforms;
    return Y;
  }
  function C(M, y) {
    let Y;
    for (let ce = 0, B = c.length; ce < B; ce++) {
      const H = c[ce];
      if (H.cacheKey === y) {
        Y = H, ++Y.usedTimes;
        break;
      }
    }
    return Y === void 0 && (Y = new zd(i, y, M, s), c.push(Y)), Y;
  }
  function L(M) {
    if (--M.usedTimes === 0) {
      const y = c.indexOf(M);
      c[y] = c[c.length - 1], c.pop(), M.destroy();
    }
  }
  function w(M) {
    l.remove(M);
  }
  function V() {
    l.dispose();
  }
  return {
    getParameters: p,
    getProgramCacheKey: d,
    getUniforms: T,
    acquireProgram: C,
    releaseProgram: L,
    releaseShaderCache: w,
    // Exposed for resource monitoring & error feedback via renderer.info:
    programs: c,
    dispose: V
  };
}
function Wd() {
  let i = /* @__PURE__ */ new WeakMap();
  function e(s) {
    let o = i.get(s);
    return o === void 0 && (o = {}, i.set(s, o)), o;
  }
  function t(s) {
    i.delete(s);
  }
  function n(s, o, a) {
    i.get(s)[o] = a;
  }
  function r() {
    i = /* @__PURE__ */ new WeakMap();
  }
  return {
    get: e,
    remove: t,
    update: n,
    dispose: r
  };
}
function Xd(i, e) {
  return i.groupOrder !== e.groupOrder ? i.groupOrder - e.groupOrder : i.renderOrder !== e.renderOrder ? i.renderOrder - e.renderOrder : i.material.id !== e.material.id ? i.material.id - e.material.id : i.z !== e.z ? i.z - e.z : i.id - e.id;
}
function Oa(i, e) {
  return i.groupOrder !== e.groupOrder ? i.groupOrder - e.groupOrder : i.renderOrder !== e.renderOrder ? i.renderOrder - e.renderOrder : i.z !== e.z ? e.z - i.z : i.id - e.id;
}
function Fa() {
  const i = [];
  let e = 0;
  const t = [], n = [], r = [];
  function s() {
    e = 0, t.length = 0, n.length = 0, r.length = 0;
  }
  function o(f, u, m, g, x, p) {
    let d = i[e];
    return d === void 0 ? (d = {
      id: f.id,
      object: f,
      geometry: u,
      material: m,
      groupOrder: g,
      renderOrder: f.renderOrder,
      z: x,
      group: p
    }, i[e] = d) : (d.id = f.id, d.object = f, d.geometry = u, d.material = m, d.groupOrder = g, d.renderOrder = f.renderOrder, d.z = x, d.group = p), e++, d;
  }
  function a(f, u, m, g, x, p) {
    const d = o(f, u, m, g, x, p);
    m.transmission > 0 ? n.push(d) : m.transparent === !0 ? r.push(d) : t.push(d);
  }
  function l(f, u, m, g, x, p) {
    const d = o(f, u, m, g, x, p);
    m.transmission > 0 ? n.unshift(d) : m.transparent === !0 ? r.unshift(d) : t.unshift(d);
  }
  function c(f, u) {
    t.length > 1 && t.sort(f || Xd), n.length > 1 && n.sort(u || Oa), r.length > 1 && r.sort(u || Oa);
  }
  function h() {
    for (let f = e, u = i.length; f < u; f++) {
      const m = i[f];
      if (m.id === null)
        break;
      m.id = null, m.object = null, m.geometry = null, m.material = null, m.group = null;
    }
  }
  return {
    opaque: t,
    transmissive: n,
    transparent: r,
    init: s,
    push: a,
    unshift: l,
    finish: h,
    sort: c
  };
}
function Yd() {
  let i = /* @__PURE__ */ new WeakMap();
  function e(n, r) {
    const s = i.get(n);
    let o;
    return s === void 0 ? (o = new Fa(), i.set(n, [o])) : r >= s.length ? (o = new Fa(), s.push(o)) : o = s[r], o;
  }
  function t() {
    i = /* @__PURE__ */ new WeakMap();
  }
  return {
    get: e,
    dispose: t
  };
}
function qd() {
  const i = {};
  return {
    get: function(e) {
      if (i[e.id] !== void 0)
        return i[e.id];
      let t;
      switch (e.type) {
        case "DirectionalLight":
          t = {
            direction: new U(),
            color: new We()
          };
          break;
        case "SpotLight":
          t = {
            position: new U(),
            direction: new U(),
            color: new We(),
            distance: 0,
            coneCos: 0,
            penumbraCos: 0,
            decay: 0
          };
          break;
        case "PointLight":
          t = {
            position: new U(),
            color: new We(),
            distance: 0,
            decay: 0
          };
          break;
        case "HemisphereLight":
          t = {
            direction: new U(),
            skyColor: new We(),
            groundColor: new We()
          };
          break;
        case "RectAreaLight":
          t = {
            color: new We(),
            position: new U(),
            halfWidth: new U(),
            halfHeight: new U()
          };
          break;
      }
      return i[e.id] = t, t;
    }
  };
}
function jd() {
  const i = {};
  return {
    get: function(e) {
      if (i[e.id] !== void 0)
        return i[e.id];
      let t;
      switch (e.type) {
        case "DirectionalLight":
          t = {
            shadowBias: 0,
            shadowNormalBias: 0,
            shadowRadius: 1,
            shadowMapSize: new oe()
          };
          break;
        case "SpotLight":
          t = {
            shadowBias: 0,
            shadowNormalBias: 0,
            shadowRadius: 1,
            shadowMapSize: new oe()
          };
          break;
        case "PointLight":
          t = {
            shadowBias: 0,
            shadowNormalBias: 0,
            shadowRadius: 1,
            shadowMapSize: new oe(),
            shadowCameraNear: 1,
            shadowCameraFar: 1e3
          };
          break;
      }
      return i[e.id] = t, t;
    }
  };
}
let Zd = 0;
function Kd(i, e) {
  return (e.castShadow ? 2 : 0) - (i.castShadow ? 2 : 0) + (e.map ? 1 : 0) - (i.map ? 1 : 0);
}
function Jd(i, e) {
  const t = new qd(), n = jd(), r = {
    version: 0,
    hash: {
      directionalLength: -1,
      pointLength: -1,
      spotLength: -1,
      rectAreaLength: -1,
      hemiLength: -1,
      numDirectionalShadows: -1,
      numPointShadows: -1,
      numSpotShadows: -1,
      numSpotMaps: -1
    },
    ambient: [0, 0, 0],
    probe: [],
    directional: [],
    directionalShadow: [],
    directionalShadowMap: [],
    directionalShadowMatrix: [],
    spot: [],
    spotLightMap: [],
    spotShadow: [],
    spotShadowMap: [],
    spotLightMatrix: [],
    rectArea: [],
    rectAreaLTC1: null,
    rectAreaLTC2: null,
    point: [],
    pointShadow: [],
    pointShadowMap: [],
    pointShadowMatrix: [],
    hemi: [],
    numSpotLightShadowsWithMaps: 0
  };
  for (let h = 0; h < 9; h++)
    r.probe.push(new U());
  const s = new U(), o = new nt(), a = new nt();
  function l(h, f) {
    let u = 0, m = 0, g = 0;
    for (let Y = 0; Y < 9; Y++)
      r.probe[Y].set(0, 0, 0);
    let x = 0, p = 0, d = 0, A = 0, _ = 0, T = 0, C = 0, L = 0, w = 0, V = 0;
    h.sort(Kd);
    const M = f === !0 ? Math.PI : 1;
    for (let Y = 0, ce = h.length; Y < ce; Y++) {
      const B = h[Y], H = B.color, G = B.intensity, Q = B.distance, X = B.shadow && B.shadow.map ? B.shadow.map.texture : null;
      if (B.isAmbientLight)
        u += H.r * G * M, m += H.g * G * M, g += H.b * G * M;
      else if (B.isLightProbe)
        for (let j = 0; j < 9; j++)
          r.probe[j].addScaledVector(B.sh.coefficients[j], G);
      else if (B.isDirectionalLight) {
        const j = t.get(B);
        if (j.color.copy(B.color).multiplyScalar(B.intensity * M), B.castShadow) {
          const J = B.shadow, ee = n.get(B);
          ee.shadowBias = J.bias, ee.shadowNormalBias = J.normalBias, ee.shadowRadius = J.radius, ee.shadowMapSize = J.mapSize, r.directionalShadow[x] = ee, r.directionalShadowMap[x] = X, r.directionalShadowMatrix[x] = B.shadow.matrix, T++;
        }
        r.directional[x] = j, x++;
      } else if (B.isSpotLight) {
        const j = t.get(B);
        j.position.setFromMatrixPosition(B.matrixWorld), j.color.copy(H).multiplyScalar(G * M), j.distance = Q, j.coneCos = Math.cos(B.angle), j.penumbraCos = Math.cos(B.angle * (1 - B.penumbra)), j.decay = B.decay, r.spot[d] = j;
        const J = B.shadow;
        if (B.map && (r.spotLightMap[w] = B.map, w++, J.updateMatrices(B), B.castShadow && V++), r.spotLightMatrix[d] = J.matrix, B.castShadow) {
          const ee = n.get(B);
          ee.shadowBias = J.bias, ee.shadowNormalBias = J.normalBias, ee.shadowRadius = J.radius, ee.shadowMapSize = J.mapSize, r.spotShadow[d] = ee, r.spotShadowMap[d] = X, L++;
        }
        d++;
      } else if (B.isRectAreaLight) {
        const j = t.get(B);
        j.color.copy(H).multiplyScalar(G), j.halfWidth.set(B.width * 0.5, 0, 0), j.halfHeight.set(0, B.height * 0.5, 0), r.rectArea[A] = j, A++;
      } else if (B.isPointLight) {
        const j = t.get(B);
        if (j.color.copy(B.color).multiplyScalar(B.intensity * M), j.distance = B.distance, j.decay = B.decay, B.castShadow) {
          const J = B.shadow, ee = n.get(B);
          ee.shadowBias = J.bias, ee.shadowNormalBias = J.normalBias, ee.shadowRadius = J.radius, ee.shadowMapSize = J.mapSize, ee.shadowCameraNear = J.camera.near, ee.shadowCameraFar = J.camera.far, r.pointShadow[p] = ee, r.pointShadowMap[p] = X, r.pointShadowMatrix[p] = B.shadow.matrix, C++;
        }
        r.point[p] = j, p++;
      } else if (B.isHemisphereLight) {
        const j = t.get(B);
        j.skyColor.copy(B.color).multiplyScalar(G * M), j.groundColor.copy(B.groundColor).multiplyScalar(G * M), r.hemi[_] = j, _++;
      }
    }
    A > 0 && (e.isWebGL2 || i.has("OES_texture_float_linear") === !0 ? (r.rectAreaLTC1 = ue.LTC_FLOAT_1, r.rectAreaLTC2 = ue.LTC_FLOAT_2) : i.has("OES_texture_half_float_linear") === !0 ? (r.rectAreaLTC1 = ue.LTC_HALF_1, r.rectAreaLTC2 = ue.LTC_HALF_2) : console.error("THREE.WebGLRenderer: Unable to use RectAreaLight. Missing WebGL extensions.")), r.ambient[0] = u, r.ambient[1] = m, r.ambient[2] = g;
    const y = r.hash;
    (y.directionalLength !== x || y.pointLength !== p || y.spotLength !== d || y.rectAreaLength !== A || y.hemiLength !== _ || y.numDirectionalShadows !== T || y.numPointShadows !== C || y.numSpotShadows !== L || y.numSpotMaps !== w) && (r.directional.length = x, r.spot.length = d, r.rectArea.length = A, r.point.length = p, r.hemi.length = _, r.directionalShadow.length = T, r.directionalShadowMap.length = T, r.pointShadow.length = C, r.pointShadowMap.length = C, r.spotShadow.length = L, r.spotShadowMap.length = L, r.directionalShadowMatrix.length = T, r.pointShadowMatrix.length = C, r.spotLightMatrix.length = L + w - V, r.spotLightMap.length = w, r.numSpotLightShadowsWithMaps = V, y.directionalLength = x, y.pointLength = p, y.spotLength = d, y.rectAreaLength = A, y.hemiLength = _, y.numDirectionalShadows = T, y.numPointShadows = C, y.numSpotShadows = L, y.numSpotMaps = w, r.version = Zd++);
  }
  function c(h, f) {
    let u = 0, m = 0, g = 0, x = 0, p = 0;
    const d = f.matrixWorldInverse;
    for (let A = 0, _ = h.length; A < _; A++) {
      const T = h[A];
      if (T.isDirectionalLight) {
        const C = r.directional[u];
        C.direction.setFromMatrixPosition(T.matrixWorld), s.setFromMatrixPosition(T.target.matrixWorld), C.direction.sub(s), C.direction.transformDirection(d), u++;
      } else if (T.isSpotLight) {
        const C = r.spot[g];
        C.position.setFromMatrixPosition(T.matrixWorld), C.position.applyMatrix4(d), C.direction.setFromMatrixPosition(T.matrixWorld), s.setFromMatrixPosition(T.target.matrixWorld), C.direction.sub(s), C.direction.transformDirection(d), g++;
      } else if (T.isRectAreaLight) {
        const C = r.rectArea[x];
        C.position.setFromMatrixPosition(T.matrixWorld), C.position.applyMatrix4(d), a.identity(), o.copy(T.matrixWorld), o.premultiply(d), a.extractRotation(o), C.halfWidth.set(T.width * 0.5, 0, 0), C.halfHeight.set(0, T.height * 0.5, 0), C.halfWidth.applyMatrix4(a), C.halfHeight.applyMatrix4(a), x++;
      } else if (T.isPointLight) {
        const C = r.point[m];
        C.position.setFromMatrixPosition(T.matrixWorld), C.position.applyMatrix4(d), m++;
      } else if (T.isHemisphereLight) {
        const C = r.hemi[p];
        C.direction.setFromMatrixPosition(T.matrixWorld), C.direction.transformDirection(d), p++;
      }
    }
  }
  return {
    setup: l,
    setupView: c,
    state: r
  };
}
function Ba(i, e) {
  const t = new Jd(i, e), n = [], r = [];
  function s() {
    n.length = 0, r.length = 0;
  }
  function o(f) {
    n.push(f);
  }
  function a(f) {
    r.push(f);
  }
  function l(f) {
    t.setup(n, f);
  }
  function c(f) {
    t.setupView(n, f);
  }
  return {
    init: s,
    state: {
      lightsArray: n,
      shadowsArray: r,
      lights: t
    },
    setupLights: l,
    setupLightsView: c,
    pushLight: o,
    pushShadow: a
  };
}
function $d(i, e) {
  let t = /* @__PURE__ */ new WeakMap();
  function n(s, o = 0) {
    const a = t.get(s);
    let l;
    return a === void 0 ? (l = new Ba(i, e), t.set(s, [l])) : o >= a.length ? (l = new Ba(i, e), a.push(l)) : l = a[o], l;
  }
  function r() {
    t = /* @__PURE__ */ new WeakMap();
  }
  return {
    get: n,
    dispose: r
  };
}
class Qd extends Ai {
  constructor(e) {
    super(), this.isMeshDepthMaterial = !0, this.type = "MeshDepthMaterial", this.depthPacking = Il, this.map = null, this.alphaMap = null, this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.wireframe = !1, this.wireframeLinewidth = 1, this.setValues(e);
  }
  copy(e) {
    return super.copy(e), this.depthPacking = e.depthPacking, this.map = e.map, this.alphaMap = e.alphaMap, this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this;
  }
}
class ep extends Ai {
  constructor(e) {
    super(), this.isMeshDistanceMaterial = !0, this.type = "MeshDistanceMaterial", this.map = null, this.alphaMap = null, this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.setValues(e);
  }
  copy(e) {
    return super.copy(e), this.map = e.map, this.alphaMap = e.alphaMap, this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this;
  }
}
const tp = `void main() {
	gl_Position = vec4( position, 1.0 );
}`, np = `uniform sampler2D shadow_pass;
uniform vec2 resolution;
uniform float radius;
#include <packing>
void main() {
	const float samples = float( VSM_SAMPLES );
	float mean = 0.0;
	float squared_mean = 0.0;
	float uvStride = samples <= 1.0 ? 0.0 : 2.0 / ( samples - 1.0 );
	float uvStart = samples <= 1.0 ? 0.0 : - 1.0;
	for ( float i = 0.0; i < samples; i ++ ) {
		float uvOffset = uvStart + i * uvStride;
		#ifdef HORIZONTAL_PASS
			vec2 distribution = unpackRGBATo2Half( texture2D( shadow_pass, ( gl_FragCoord.xy + vec2( uvOffset, 0.0 ) * radius ) / resolution ) );
			mean += distribution.x;
			squared_mean += distribution.y * distribution.y + distribution.x * distribution.x;
		#else
			float depth = unpackRGBAToDepth( texture2D( shadow_pass, ( gl_FragCoord.xy + vec2( 0.0, uvOffset ) * radius ) / resolution ) );
			mean += depth;
			squared_mean += depth * depth;
		#endif
	}
	mean = mean / samples;
	squared_mean = squared_mean / samples;
	float std_dev = sqrt( squared_mean - mean * mean );
	gl_FragColor = pack2HalfToRGBA( vec2( mean, std_dev ) );
}`;
function ip(i, e, t) {
  let n = new us();
  const r = new oe(), s = new oe(), o = new ot(), a = new Qd({ depthPacking: Nl }), l = new ep(), c = {}, h = t.maxTextureSize, f = { [ln]: gt, [gt]: ln, [Mt]: Mt }, u = new bn({
    defines: {
      VSM_SAMPLES: 8
    },
    uniforms: {
      shadow_pass: { value: null },
      resolution: { value: new oe() },
      radius: { value: 4 }
    },
    vertexShader: tp,
    fragmentShader: np
  }), m = u.clone();
  m.defines.HORIZONTAL_PASS = 1;
  const g = new cn();
  g.setAttribute(
    "position",
    new Ot(
      new Float32Array([-1, -1, 0.5, 3, -1, 0.5, -1, 3, 0.5]),
      3
    )
  );
  const x = new Nt(g, u), p = this;
  this.enabled = !1, this.autoUpdate = !0, this.needsUpdate = !1, this.type = $a;
  let d = this.type;
  this.render = function(C, L, w) {
    if (p.enabled === !1 || p.autoUpdate === !1 && p.needsUpdate === !1 || C.length === 0)
      return;
    const V = i.getRenderTarget(), M = i.getActiveCubeFace(), y = i.getActiveMipmapLevel(), Y = i.state;
    Y.setBlending(sn), Y.buffers.color.setClear(1, 1, 1, 1), Y.buffers.depth.setTest(!0), Y.setScissorTest(!1);
    const ce = d !== Yt && this.type === Yt, B = d === Yt && this.type !== Yt;
    for (let H = 0, G = C.length; H < G; H++) {
      const Q = C[H], X = Q.shadow;
      if (X === void 0) {
        console.warn("THREE.WebGLShadowMap:", Q, "has no shadow.");
        continue;
      }
      if (X.autoUpdate === !1 && X.needsUpdate === !1)
        continue;
      r.copy(X.mapSize);
      const j = X.getFrameExtents();
      if (r.multiply(j), s.copy(X.mapSize), (r.x > h || r.y > h) && (r.x > h && (s.x = Math.floor(h / j.x), r.x = s.x * j.x, X.mapSize.x = s.x), r.y > h && (s.y = Math.floor(h / j.y), r.y = s.y * j.y, X.mapSize.y = s.y)), X.map === null || ce === !0 || B === !0) {
        const ee = this.type !== Yt ? { minFilter: mt, magFilter: mt } : {};
        X.map !== null && X.map.dispose(), X.map = new Tn(r.x, r.y, ee), X.map.texture.name = Q.name + ".shadowMap", X.camera.updateProjectionMatrix();
      }
      i.setRenderTarget(X.map), i.clear();
      const J = X.getViewportCount();
      for (let ee = 0; ee < J; ee++) {
        const I = X.getViewport(ee);
        o.set(
          s.x * I.x,
          s.y * I.y,
          s.x * I.z,
          s.y * I.w
        ), Y.viewport(o), X.updateMatrices(Q, ee), n = X.getFrustum(), T(L, w, X.camera, Q, this.type);
      }
      X.isPointLightShadow !== !0 && this.type === Yt && A(X, w), X.needsUpdate = !1;
    }
    d = this.type, p.needsUpdate = !1, i.setRenderTarget(V, M, y);
  };
  function A(C, L) {
    const w = e.update(x);
    u.defines.VSM_SAMPLES !== C.blurSamples && (u.defines.VSM_SAMPLES = C.blurSamples, m.defines.VSM_SAMPLES = C.blurSamples, u.needsUpdate = !0, m.needsUpdate = !0), C.mapPass === null && (C.mapPass = new Tn(r.x, r.y)), u.uniforms.shadow_pass.value = C.map.texture, u.uniforms.resolution.value = C.mapSize, u.uniforms.radius.value = C.radius, i.setRenderTarget(C.mapPass), i.clear(), i.renderBufferDirect(L, null, w, u, x, null), m.uniforms.shadow_pass.value = C.mapPass.texture, m.uniforms.resolution.value = C.mapSize, m.uniforms.radius.value = C.radius, i.setRenderTarget(C.map), i.clear(), i.renderBufferDirect(L, null, w, m, x, null);
  }
  function _(C, L, w, V) {
    let M = null;
    const y = w.isPointLight === !0 ? C.customDistanceMaterial : C.customDepthMaterial;
    if (y !== void 0)
      M = y;
    else if (M = w.isPointLight === !0 ? l : a, i.localClippingEnabled && L.clipShadows === !0 && Array.isArray(L.clippingPlanes) && L.clippingPlanes.length !== 0 || L.displacementMap && L.displacementScale !== 0 || L.alphaMap && L.alphaTest > 0 || L.map && L.alphaTest > 0) {
      const Y = M.uuid, ce = L.uuid;
      let B = c[Y];
      B === void 0 && (B = {}, c[Y] = B);
      let H = B[ce];
      H === void 0 && (H = M.clone(), B[ce] = H), M = H;
    }
    if (M.visible = L.visible, M.wireframe = L.wireframe, V === Yt ? M.side = L.shadowSide !== null ? L.shadowSide : L.side : M.side = L.shadowSide !== null ? L.shadowSide : f[L.side], M.alphaMap = L.alphaMap, M.alphaTest = L.alphaTest, M.map = L.map, M.clipShadows = L.clipShadows, M.clippingPlanes = L.clippingPlanes, M.clipIntersection = L.clipIntersection, M.displacementMap = L.displacementMap, M.displacementScale = L.displacementScale, M.displacementBias = L.displacementBias, M.wireframeLinewidth = L.wireframeLinewidth, M.linewidth = L.linewidth, w.isPointLight === !0 && M.isMeshDistanceMaterial === !0) {
      const Y = i.properties.get(M);
      Y.light = w;
    }
    return M;
  }
  function T(C, L, w, V, M) {
    if (C.visible === !1)
      return;
    if (C.layers.test(L.layers) && (C.isMesh || C.isLine || C.isPoints) && (C.castShadow || C.receiveShadow && M === Yt) && (!C.frustumCulled || n.intersectsObject(C))) {
      C.modelViewMatrix.multiplyMatrices(w.matrixWorldInverse, C.matrixWorld);
      const ce = e.update(C), B = C.material;
      if (Array.isArray(B)) {
        const H = ce.groups;
        for (let G = 0, Q = H.length; G < Q; G++) {
          const X = H[G], j = B[X.materialIndex];
          if (j && j.visible) {
            const J = _(C, j, V, M);
            i.renderBufferDirect(w, null, ce, J, C, X);
          }
        }
      } else if (B.visible) {
        const H = _(C, B, V, M);
        i.renderBufferDirect(w, null, ce, H, C, null);
      }
    }
    const Y = C.children;
    for (let ce = 0, B = Y.length; ce < B; ce++)
      T(Y[ce], L, w, V, M);
  }
}
function rp(i, e, t) {
  const n = t.isWebGL2;
  function r() {
    let P = !1;
    const fe = new ot();
    let F = null;
    const ne = new ot(0, 0, 0, 0);
    return {
      setMask: function(de) {
        F !== de && !P && (i.colorMask(de, de, de, de), F = de);
      },
      setLocked: function(de) {
        P = de;
      },
      setClear: function(de, Fe, Xe, qe, Zt) {
        Zt === !0 && (de *= qe, Fe *= qe, Xe *= qe), fe.set(de, Fe, Xe, qe), ne.equals(fe) === !1 && (i.clearColor(de, Fe, Xe, qe), ne.copy(fe));
      },
      reset: function() {
        P = !1, F = null, ne.set(-1, 0, 0, 0);
      }
    };
  }
  function s() {
    let P = !1, fe = null, F = null, ne = null;
    return {
      setTest: function(de) {
        de ? be(i.DEPTH_TEST) : Te(i.DEPTH_TEST);
      },
      setMask: function(de) {
        fe !== de && !P && (i.depthMask(de), fe = de);
      },
      setFunc: function(de) {
        if (F !== de) {
          switch (de) {
            case hl:
              i.depthFunc(i.NEVER);
              break;
            case ul:
              i.depthFunc(i.ALWAYS);
              break;
            case fl:
              i.depthFunc(i.LESS);
              break;
            case Yr:
              i.depthFunc(i.LEQUAL);
              break;
            case dl:
              i.depthFunc(i.EQUAL);
              break;
            case pl:
              i.depthFunc(i.GEQUAL);
              break;
            case ml:
              i.depthFunc(i.GREATER);
              break;
            case gl:
              i.depthFunc(i.NOTEQUAL);
              break;
            default:
              i.depthFunc(i.LEQUAL);
          }
          F = de;
        }
      },
      setLocked: function(de) {
        P = de;
      },
      setClear: function(de) {
        ne !== de && (i.clearDepth(de), ne = de);
      },
      reset: function() {
        P = !1, fe = null, F = null, ne = null;
      }
    };
  }
  function o() {
    let P = !1, fe = null, F = null, ne = null, de = null, Fe = null, Xe = null, qe = null, Zt = null;
    return {
      setTest: function(je) {
        P || (je ? be(i.STENCIL_TEST) : Te(i.STENCIL_TEST));
      },
      setMask: function(je) {
        fe !== je && !P && (i.stencilMask(je), fe = je);
      },
      setFunc: function(je, Dt, ut) {
        (F !== je || ne !== Dt || de !== ut) && (i.stencilFunc(je, Dt, ut), F = je, ne = Dt, de = ut);
      },
      setOp: function(je, Dt, ut) {
        (Fe !== je || Xe !== Dt || qe !== ut) && (i.stencilOp(je, Dt, ut), Fe = je, Xe = Dt, qe = ut);
      },
      setLocked: function(je) {
        P = je;
      },
      setClear: function(je) {
        Zt !== je && (i.clearStencil(je), Zt = je);
      },
      reset: function() {
        P = !1, fe = null, F = null, ne = null, de = null, Fe = null, Xe = null, qe = null, Zt = null;
      }
    };
  }
  const a = new r(), l = new s(), c = new o(), h = /* @__PURE__ */ new WeakMap(), f = /* @__PURE__ */ new WeakMap();
  let u = {}, m = {}, g = /* @__PURE__ */ new WeakMap(), x = [], p = null, d = !1, A = null, _ = null, T = null, C = null, L = null, w = null, V = null, M = !1, y = null, Y = null, ce = null, B = null, H = null;
  const G = i.getParameter(i.MAX_COMBINED_TEXTURE_IMAGE_UNITS);
  let Q = !1, X = 0;
  const j = i.getParameter(i.VERSION);
  j.indexOf("WebGL") !== -1 ? (X = parseFloat(/^WebGL (\d)/.exec(j)[1]), Q = X >= 1) : j.indexOf("OpenGL ES") !== -1 && (X = parseFloat(/^OpenGL ES (\d)/.exec(j)[1]), Q = X >= 2);
  let J = null, ee = {};
  const I = i.getParameter(i.SCISSOR_BOX), q = i.getParameter(i.VIEWPORT), pe = new ot().fromArray(I), me = new ot().fromArray(q);
  function ve(P, fe, F, ne) {
    const de = new Uint8Array(4), Fe = i.createTexture();
    i.bindTexture(P, Fe), i.texParameteri(P, i.TEXTURE_MIN_FILTER, i.NEAREST), i.texParameteri(P, i.TEXTURE_MAG_FILTER, i.NEAREST);
    for (let Xe = 0; Xe < F; Xe++)
      n && (P === i.TEXTURE_3D || P === i.TEXTURE_2D_ARRAY) ? i.texImage3D(fe, 0, i.RGBA, 1, 1, ne, 0, i.RGBA, i.UNSIGNED_BYTE, de) : i.texImage2D(fe + Xe, 0, i.RGBA, 1, 1, 0, i.RGBA, i.UNSIGNED_BYTE, de);
    return Fe;
  }
  const Ae = {};
  Ae[i.TEXTURE_2D] = ve(i.TEXTURE_2D, i.TEXTURE_2D, 1), Ae[i.TEXTURE_CUBE_MAP] = ve(i.TEXTURE_CUBE_MAP, i.TEXTURE_CUBE_MAP_POSITIVE_X, 6), n && (Ae[i.TEXTURE_2D_ARRAY] = ve(i.TEXTURE_2D_ARRAY, i.TEXTURE_2D_ARRAY, 1, 1), Ae[i.TEXTURE_3D] = ve(i.TEXTURE_3D, i.TEXTURE_3D, 1, 1)), a.setClear(0, 0, 0, 1), l.setClear(1), c.setClear(0), be(i.DEPTH_TEST), l.setFunc(Yr), $(!1), ye(ys), be(i.CULL_FACE), Z(sn);
  function be(P) {
    u[P] !== !0 && (i.enable(P), u[P] = !0);
  }
  function Te(P) {
    u[P] !== !1 && (i.disable(P), u[P] = !1);
  }
  function ke(P, fe) {
    return m[P] !== fe ? (i.bindFramebuffer(P, fe), m[P] = fe, n && (P === i.DRAW_FRAMEBUFFER && (m[i.FRAMEBUFFER] = fe), P === i.FRAMEBUFFER && (m[i.DRAW_FRAMEBUFFER] = fe)), !0) : !1;
  }
  function Ye(P, fe) {
    let F = x, ne = !1;
    if (P)
      if (F = g.get(fe), F === void 0 && (F = [], g.set(fe, F)), P.isWebGLMultipleRenderTargets) {
        const de = P.texture;
        if (F.length !== de.length || F[0] !== i.COLOR_ATTACHMENT0) {
          for (let Fe = 0, Xe = de.length; Fe < Xe; Fe++)
            F[Fe] = i.COLOR_ATTACHMENT0 + Fe;
          F.length = de.length, ne = !0;
        }
      } else
        F[0] !== i.COLOR_ATTACHMENT0 && (F[0] = i.COLOR_ATTACHMENT0, ne = !0);
    else
      F[0] !== i.BACK && (F[0] = i.BACK, ne = !0);
    ne && (t.isWebGL2 ? i.drawBuffers(F) : e.get("WEBGL_draw_buffers").drawBuffersWEBGL(F));
  }
  function we(P) {
    return p !== P ? (i.useProgram(P), p = P, !0) : !1;
  }
  const b = {
    [jn]: i.FUNC_ADD,
    [Qo]: i.FUNC_SUBTRACT,
    [el]: i.FUNC_REVERSE_SUBTRACT
  };
  if (n)
    b[ws] = i.MIN, b[Rs] = i.MAX;
  else {
    const P = e.get("EXT_blend_minmax");
    P !== null && (b[ws] = P.MIN_EXT, b[Rs] = P.MAX_EXT);
  }
  const le = {
    [tl]: i.ZERO,
    [nl]: i.ONE,
    [il]: i.SRC_COLOR,
    [Qa]: i.SRC_ALPHA,
    [cl]: i.SRC_ALPHA_SATURATE,
    [ol]: i.DST_COLOR,
    [sl]: i.DST_ALPHA,
    [rl]: i.ONE_MINUS_SRC_COLOR,
    [eo]: i.ONE_MINUS_SRC_ALPHA,
    [ll]: i.ONE_MINUS_DST_COLOR,
    [al]: i.ONE_MINUS_DST_ALPHA
  };
  function Z(P, fe, F, ne, de, Fe, Xe, qe) {
    if (P === sn) {
      d === !0 && (Te(i.BLEND), d = !1);
      return;
    }
    if (d === !1 && (be(i.BLEND), d = !0), P !== $o) {
      if (P !== A || qe !== M) {
        if ((_ !== jn || L !== jn) && (i.blendEquation(i.FUNC_ADD), _ = jn, L = jn), qe)
          switch (P) {
            case $n:
              i.blendFuncSeparate(i.ONE, i.ONE_MINUS_SRC_ALPHA, i.ONE, i.ONE_MINUS_SRC_ALPHA);
              break;
            case Ts:
              i.blendFunc(i.ONE, i.ONE);
              break;
            case As:
              i.blendFuncSeparate(i.ZERO, i.ONE_MINUS_SRC_COLOR, i.ZERO, i.ONE);
              break;
            case bs:
              i.blendFuncSeparate(i.ZERO, i.SRC_COLOR, i.ZERO, i.SRC_ALPHA);
              break;
            default:
              console.error("THREE.WebGLState: Invalid blending: ", P);
              break;
          }
        else
          switch (P) {
            case $n:
              i.blendFuncSeparate(i.SRC_ALPHA, i.ONE_MINUS_SRC_ALPHA, i.ONE, i.ONE_MINUS_SRC_ALPHA);
              break;
            case Ts:
              i.blendFunc(i.SRC_ALPHA, i.ONE);
              break;
            case As:
              i.blendFuncSeparate(i.ZERO, i.ONE_MINUS_SRC_COLOR, i.ZERO, i.ONE);
              break;
            case bs:
              i.blendFunc(i.ZERO, i.SRC_COLOR);
              break;
            default:
              console.error("THREE.WebGLState: Invalid blending: ", P);
              break;
          }
        T = null, C = null, w = null, V = null, A = P, M = qe;
      }
      return;
    }
    de = de || fe, Fe = Fe || F, Xe = Xe || ne, (fe !== _ || de !== L) && (i.blendEquationSeparate(b[fe], b[de]), _ = fe, L = de), (F !== T || ne !== C || Fe !== w || Xe !== V) && (i.blendFuncSeparate(le[F], le[ne], le[Fe], le[Xe]), T = F, C = ne, w = Fe, V = Xe), A = P, M = !1;
  }
  function re(P, fe) {
    P.side === Mt ? Te(i.CULL_FACE) : be(i.CULL_FACE);
    let F = P.side === gt;
    fe && (F = !F), $(F), P.blending === $n && P.transparent === !1 ? Z(sn) : Z(P.blending, P.blendEquation, P.blendSrc, P.blendDst, P.blendEquationAlpha, P.blendSrcAlpha, P.blendDstAlpha, P.premultipliedAlpha), l.setFunc(P.depthFunc), l.setTest(P.depthTest), l.setMask(P.depthWrite), a.setMask(P.colorWrite);
    const ne = P.stencilWrite;
    c.setTest(ne), ne && (c.setMask(P.stencilWriteMask), c.setFunc(P.stencilFunc, P.stencilRef, P.stencilFuncMask), c.setOp(P.stencilFail, P.stencilZFail, P.stencilZPass)), Me(P.polygonOffset, P.polygonOffsetFactor, P.polygonOffsetUnits), P.alphaToCoverage === !0 ? be(i.SAMPLE_ALPHA_TO_COVERAGE) : Te(i.SAMPLE_ALPHA_TO_COVERAGE);
  }
  function $(P) {
    y !== P && (P ? i.frontFace(i.CW) : i.frontFace(i.CCW), y = P);
  }
  function ye(P) {
    P !== Zo ? (be(i.CULL_FACE), P !== Y && (P === ys ? i.cullFace(i.BACK) : P === Ko ? i.cullFace(i.FRONT) : i.cullFace(i.FRONT_AND_BACK))) : Te(i.CULL_FACE), Y = P;
  }
  function xe(P) {
    P !== ce && (Q && i.lineWidth(P), ce = P);
  }
  function Me(P, fe, F) {
    P ? (be(i.POLYGON_OFFSET_FILL), (B !== fe || H !== F) && (i.polygonOffset(fe, F), B = fe, H = F)) : Te(i.POLYGON_OFFSET_FILL);
  }
  function Ce(P) {
    P ? be(i.SCISSOR_TEST) : Te(i.SCISSOR_TEST);
  }
  function ze(P) {
    P === void 0 && (P = i.TEXTURE0 + G - 1), J !== P && (i.activeTexture(P), J = P);
  }
  function Ze(P, fe, F) {
    F === void 0 && (J === null ? F = i.TEXTURE0 + G - 1 : F = J);
    let ne = ee[F];
    ne === void 0 && (ne = { type: void 0, texture: void 0 }, ee[F] = ne), (ne.type !== P || ne.texture !== fe) && (J !== F && (i.activeTexture(F), J = F), i.bindTexture(P, fe || Ae[P]), ne.type = P, ne.texture = fe);
  }
  function E() {
    const P = ee[J];
    P !== void 0 && P.type !== void 0 && (i.bindTexture(P.type, null), P.type = void 0, P.texture = void 0);
  }
  function v() {
    try {
      i.compressedTexImage2D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function O() {
    try {
      i.compressedTexImage3D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function ie() {
    try {
      i.texSubImage2D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function te() {
    try {
      i.texSubImage3D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function se() {
    try {
      i.compressedTexSubImage2D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function Ee() {
    try {
      i.compressedTexSubImage3D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function ae() {
    try {
      i.texStorage2D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function z() {
    try {
      i.texStorage3D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function R() {
    try {
      i.texImage2D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function K() {
    try {
      i.texImage3D.apply(i, arguments);
    } catch (P) {
      console.error("THREE.WebGLState:", P);
    }
  }
  function _e(P) {
    pe.equals(P) === !1 && (i.scissor(P.x, P.y, P.z, P.w), pe.copy(P));
  }
  function he(P) {
    me.equals(P) === !1 && (i.viewport(P.x, P.y, P.z, P.w), me.copy(P));
  }
  function ge(P, fe) {
    let F = f.get(fe);
    F === void 0 && (F = /* @__PURE__ */ new WeakMap(), f.set(fe, F));
    let ne = F.get(P);
    ne === void 0 && (ne = i.getUniformBlockIndex(fe, P.name), F.set(P, ne));
  }
  function De(P, fe) {
    const ne = f.get(fe).get(P);
    h.get(fe) !== ne && (i.uniformBlockBinding(fe, ne, P.__bindingPointIndex), h.set(fe, ne));
  }
  function Ve() {
    i.disable(i.BLEND), i.disable(i.CULL_FACE), i.disable(i.DEPTH_TEST), i.disable(i.POLYGON_OFFSET_FILL), i.disable(i.SCISSOR_TEST), i.disable(i.STENCIL_TEST), i.disable(i.SAMPLE_ALPHA_TO_COVERAGE), i.blendEquation(i.FUNC_ADD), i.blendFunc(i.ONE, i.ZERO), i.blendFuncSeparate(i.ONE, i.ZERO, i.ONE, i.ZERO), i.colorMask(!0, !0, !0, !0), i.clearColor(0, 0, 0, 0), i.depthMask(!0), i.depthFunc(i.LESS), i.clearDepth(1), i.stencilMask(4294967295), i.stencilFunc(i.ALWAYS, 0, 4294967295), i.stencilOp(i.KEEP, i.KEEP, i.KEEP), i.clearStencil(0), i.cullFace(i.BACK), i.frontFace(i.CCW), i.polygonOffset(0, 0), i.activeTexture(i.TEXTURE0), i.bindFramebuffer(i.FRAMEBUFFER, null), n === !0 && (i.bindFramebuffer(i.DRAW_FRAMEBUFFER, null), i.bindFramebuffer(i.READ_FRAMEBUFFER, null)), i.useProgram(null), i.lineWidth(1), i.scissor(0, 0, i.canvas.width, i.canvas.height), i.viewport(0, 0, i.canvas.width, i.canvas.height), u = {}, J = null, ee = {}, m = {}, g = /* @__PURE__ */ new WeakMap(), x = [], p = null, d = !1, A = null, _ = null, T = null, C = null, L = null, w = null, V = null, M = !1, y = null, Y = null, ce = null, B = null, H = null, pe.set(0, 0, i.canvas.width, i.canvas.height), me.set(0, 0, i.canvas.width, i.canvas.height), a.reset(), l.reset(), c.reset();
  }
  return {
    buffers: {
      color: a,
      depth: l,
      stencil: c
    },
    enable: be,
    disable: Te,
    bindFramebuffer: ke,
    drawBuffers: Ye,
    useProgram: we,
    setBlending: Z,
    setMaterial: re,
    setFlipSided: $,
    setCullFace: ye,
    setLineWidth: xe,
    setPolygonOffset: Me,
    setScissorTest: Ce,
    activeTexture: ze,
    bindTexture: Ze,
    unbindTexture: E,
    compressedTexImage2D: v,
    compressedTexImage3D: O,
    texImage2D: R,
    texImage3D: K,
    updateUBOMapping: ge,
    uniformBlockBinding: De,
    texStorage2D: ae,
    texStorage3D: z,
    texSubImage2D: ie,
    texSubImage3D: te,
    compressedTexSubImage2D: se,
    compressedTexSubImage3D: Ee,
    scissor: _e,
    viewport: he,
    reset: Ve
  };
}
function sp(i, e, t, n, r, s, o) {
  const a = r.isWebGL2, l = r.maxTextures, c = r.maxCubemapSize, h = r.maxTextureSize, f = r.maxSamples, u = e.has("WEBGL_multisampled_render_to_texture") ? e.get("WEBGL_multisampled_render_to_texture") : null, m = typeof navigator > "u" ? !1 : /OculusBrowser/g.test(navigator.userAgent), g = /* @__PURE__ */ new WeakMap();
  let x;
  const p = /* @__PURE__ */ new WeakMap();
  let d = !1;
  try {
    d = typeof OffscreenCanvas < "u" && new OffscreenCanvas(1, 1).getContext("2d") !== null;
  } catch {
  }
  function A(E, v) {
    return d ? (
      // eslint-disable-next-line compat/compat
      new OffscreenCanvas(E, v)
    ) : ir("canvas");
  }
  function _(E, v, O, ie) {
    let te = 1;
    if ((E.width > ie || E.height > ie) && (te = ie / Math.max(E.width, E.height)), te < 1 || v === !0)
      if (typeof HTMLImageElement < "u" && E instanceof HTMLImageElement || typeof HTMLCanvasElement < "u" && E instanceof HTMLCanvasElement || typeof ImageBitmap < "u" && E instanceof ImageBitmap) {
        const se = v ? nr : Math.floor, Ee = se(te * E.width), ae = se(te * E.height);
        x === void 0 && (x = A(Ee, ae));
        const z = O ? A(Ee, ae) : x;
        return z.width = Ee, z.height = ae, z.getContext("2d").drawImage(E, 0, 0, Ee, ae), console.warn("THREE.WebGLRenderer: Texture has been resized from (" + E.width + "x" + E.height + ") to (" + Ee + "x" + ae + ")."), z;
      } else
        return "data" in E && console.warn("THREE.WebGLRenderer: Image in DataTexture is too big (" + E.width + "x" + E.height + ")."), E;
    return E;
  }
  function T(E) {
    return $r(E.width) && $r(E.height);
  }
  function C(E) {
    return a ? !1 : E.wrapS !== Lt || E.wrapT !== Lt || E.minFilter !== mt && E.minFilter !== Tt;
  }
  function L(E, v) {
    return E.generateMipmaps && v && E.minFilter !== mt && E.minFilter !== Tt;
  }
  function w(E) {
    i.generateMipmap(E);
  }
  function V(E, v, O, ie, te = !1) {
    if (a === !1)
      return v;
    if (E !== null) {
      if (i[E] !== void 0)
        return i[E];
      console.warn("THREE.WebGLRenderer: Attempt to use non-existing WebGL internal format '" + E + "'");
    }
    let se = v;
    return v === i.RED && (O === i.FLOAT && (se = i.R32F), O === i.HALF_FLOAT && (se = i.R16F), O === i.UNSIGNED_BYTE && (se = i.R8)), v === i.RED_INTEGER && (O === i.UNSIGNED_BYTE && (se = i.R8UI), O === i.UNSIGNED_SHORT && (se = i.R16UI), O === i.UNSIGNED_INT && (se = i.R32UI), O === i.BYTE && (se = i.R8I), O === i.SHORT && (se = i.R16I), O === i.INT && (se = i.R32I)), v === i.RG && (O === i.FLOAT && (se = i.RG32F), O === i.HALF_FLOAT && (se = i.RG16F), O === i.UNSIGNED_BYTE && (se = i.RG8)), v === i.RGBA && (O === i.FLOAT && (se = i.RGBA32F), O === i.HALF_FLOAT && (se = i.RGBA16F), O === i.UNSIGNED_BYTE && (se = ie === Oe && te === !1 ? i.SRGB8_ALPHA8 : i.RGBA8), O === i.UNSIGNED_SHORT_4_4_4_4 && (se = i.RGBA4), O === i.UNSIGNED_SHORT_5_5_5_1 && (se = i.RGB5_A1)), (se === i.R16F || se === i.R32F || se === i.RG16F || se === i.RG32F || se === i.RGBA16F || se === i.RGBA32F) && e.get("EXT_color_buffer_float"), se;
  }
  function M(E, v, O) {
    return L(E, O) === !0 || E.isFramebufferTexture && E.minFilter !== mt && E.minFilter !== Tt ? Math.log2(Math.max(v.width, v.height)) + 1 : E.mipmaps !== void 0 && E.mipmaps.length > 0 ? E.mipmaps.length : E.isCompressedTexture && Array.isArray(E.image) ? v.mipmaps.length : 1;
  }
  function y(E) {
    return E === mt || E === Cs || E === dr ? i.NEAREST : i.LINEAR;
  }
  function Y(E) {
    const v = E.target;
    v.removeEventListener("dispose", Y), B(v), v.isVideoTexture && g.delete(v);
  }
  function ce(E) {
    const v = E.target;
    v.removeEventListener("dispose", ce), G(v);
  }
  function B(E) {
    const v = n.get(E);
    if (v.__webglInit === void 0)
      return;
    const O = E.source, ie = p.get(O);
    if (ie) {
      const te = ie[v.__cacheKey];
      te.usedTimes--, te.usedTimes === 0 && H(E), Object.keys(ie).length === 0 && p.delete(O);
    }
    n.remove(E);
  }
  function H(E) {
    const v = n.get(E);
    i.deleteTexture(v.__webglTexture);
    const O = E.source, ie = p.get(O);
    delete ie[v.__cacheKey], o.memory.textures--;
  }
  function G(E) {
    const v = E.texture, O = n.get(E), ie = n.get(v);
    if (ie.__webglTexture !== void 0 && (i.deleteTexture(ie.__webglTexture), o.memory.textures--), E.depthTexture && E.depthTexture.dispose(), E.isWebGLCubeRenderTarget)
      for (let te = 0; te < 6; te++) {
        if (Array.isArray(O.__webglFramebuffer[te]))
          for (let se = 0; se < O.__webglFramebuffer[te].length; se++)
            i.deleteFramebuffer(O.__webglFramebuffer[te][se]);
        else
          i.deleteFramebuffer(O.__webglFramebuffer[te]);
        O.__webglDepthbuffer && i.deleteRenderbuffer(O.__webglDepthbuffer[te]);
      }
    else {
      if (Array.isArray(O.__webglFramebuffer))
        for (let te = 0; te < O.__webglFramebuffer.length; te++)
          i.deleteFramebuffer(O.__webglFramebuffer[te]);
      else
        i.deleteFramebuffer(O.__webglFramebuffer);
      if (O.__webglDepthbuffer && i.deleteRenderbuffer(O.__webglDepthbuffer), O.__webglMultisampledFramebuffer && i.deleteFramebuffer(O.__webglMultisampledFramebuffer), O.__webglColorRenderbuffer)
        for (let te = 0; te < O.__webglColorRenderbuffer.length; te++)
          O.__webglColorRenderbuffer[te] && i.deleteRenderbuffer(O.__webglColorRenderbuffer[te]);
      O.__webglDepthRenderbuffer && i.deleteRenderbuffer(O.__webglDepthRenderbuffer);
    }
    if (E.isWebGLMultipleRenderTargets)
      for (let te = 0, se = v.length; te < se; te++) {
        const Ee = n.get(v[te]);
        Ee.__webglTexture && (i.deleteTexture(Ee.__webglTexture), o.memory.textures--), n.remove(v[te]);
      }
    n.remove(v), n.remove(E);
  }
  let Q = 0;
  function X() {
    Q = 0;
  }
  function j() {
    const E = Q;
    return E >= l && console.warn("THREE.WebGLTextures: Trying to use " + E + " texture units while this GPU supports only " + l), Q += 1, E;
  }
  function J(E) {
    const v = [];
    return v.push(E.wrapS), v.push(E.wrapT), v.push(E.wrapR || 0), v.push(E.magFilter), v.push(E.minFilter), v.push(E.anisotropy), v.push(E.internalFormat), v.push(E.format), v.push(E.type), v.push(E.generateMipmaps), v.push(E.premultiplyAlpha), v.push(E.flipY), v.push(E.unpackAlignment), v.push(E.colorSpace), v.join();
  }
  function ee(E, v) {
    const O = n.get(E);
    if (E.isVideoTexture && ze(E), E.isRenderTargetTexture === !1 && E.version > 0 && O.__version !== E.version) {
      const ie = E.image;
      if (ie === null)
        console.warn("THREE.WebGLRenderer: Texture marked for update but no image data found.");
      else if (ie.complete === !1)
        console.warn("THREE.WebGLRenderer: Texture marked for update but image is incomplete");
      else {
        ke(O, E, v);
        return;
      }
    }
    t.bindTexture(i.TEXTURE_2D, O.__webglTexture, i.TEXTURE0 + v);
  }
  function I(E, v) {
    const O = n.get(E);
    if (E.version > 0 && O.__version !== E.version) {
      ke(O, E, v);
      return;
    }
    t.bindTexture(i.TEXTURE_2D_ARRAY, O.__webglTexture, i.TEXTURE0 + v);
  }
  function q(E, v) {
    const O = n.get(E);
    if (E.version > 0 && O.__version !== E.version) {
      ke(O, E, v);
      return;
    }
    t.bindTexture(i.TEXTURE_3D, O.__webglTexture, i.TEXTURE0 + v);
  }
  function pe(E, v) {
    const O = n.get(E);
    if (E.version > 0 && O.__version !== E.version) {
      Ye(O, E, v);
      return;
    }
    t.bindTexture(i.TEXTURE_CUBE_MAP, O.__webglTexture, i.TEXTURE0 + v);
  }
  const me = {
    [Zr]: i.REPEAT,
    [Lt]: i.CLAMP_TO_EDGE,
    [Kr]: i.MIRRORED_REPEAT
  }, ve = {
    [mt]: i.NEAREST,
    [Cs]: i.NEAREST_MIPMAP_NEAREST,
    [dr]: i.NEAREST_MIPMAP_LINEAR,
    [Tt]: i.LINEAR,
    [Tl]: i.LINEAR_MIPMAP_NEAREST,
    [vi]: i.LINEAR_MIPMAP_LINEAR
  }, Ae = {
    [Bl]: i.NEVER,
    [Xl]: i.ALWAYS,
    [zl]: i.LESS,
    [Gl]: i.LEQUAL,
    [Hl]: i.EQUAL,
    [Wl]: i.GEQUAL,
    [Vl]: i.GREATER,
    [kl]: i.NOTEQUAL
  };
  function be(E, v, O) {
    if (O ? (i.texParameteri(E, i.TEXTURE_WRAP_S, me[v.wrapS]), i.texParameteri(E, i.TEXTURE_WRAP_T, me[v.wrapT]), (E === i.TEXTURE_3D || E === i.TEXTURE_2D_ARRAY) && i.texParameteri(E, i.TEXTURE_WRAP_R, me[v.wrapR]), i.texParameteri(E, i.TEXTURE_MAG_FILTER, ve[v.magFilter]), i.texParameteri(E, i.TEXTURE_MIN_FILTER, ve[v.minFilter])) : (i.texParameteri(E, i.TEXTURE_WRAP_S, i.CLAMP_TO_EDGE), i.texParameteri(E, i.TEXTURE_WRAP_T, i.CLAMP_TO_EDGE), (E === i.TEXTURE_3D || E === i.TEXTURE_2D_ARRAY) && i.texParameteri(E, i.TEXTURE_WRAP_R, i.CLAMP_TO_EDGE), (v.wrapS !== Lt || v.wrapT !== Lt) && console.warn("THREE.WebGLRenderer: Texture is not power of two. Texture.wrapS and Texture.wrapT should be set to THREE.ClampToEdgeWrapping."), i.texParameteri(E, i.TEXTURE_MAG_FILTER, y(v.magFilter)), i.texParameteri(E, i.TEXTURE_MIN_FILTER, y(v.minFilter)), v.minFilter !== mt && v.minFilter !== Tt && console.warn("THREE.WebGLRenderer: Texture is not power of two. Texture.minFilter should be set to THREE.NearestFilter or THREE.LinearFilter.")), v.compareFunction && (i.texParameteri(E, i.TEXTURE_COMPARE_MODE, i.COMPARE_REF_TO_TEXTURE), i.texParameteri(E, i.TEXTURE_COMPARE_FUNC, Ae[v.compareFunction])), e.has("EXT_texture_filter_anisotropic") === !0) {
      const ie = e.get("EXT_texture_filter_anisotropic");
      if (v.magFilter === mt || v.minFilter !== dr && v.minFilter !== vi || v.type === rn && e.has("OES_texture_float_linear") === !1 || a === !1 && v.type === xi && e.has("OES_texture_half_float_linear") === !1)
        return;
      (v.anisotropy > 1 || n.get(v).__currentAnisotropy) && (i.texParameterf(E, ie.TEXTURE_MAX_ANISOTROPY_EXT, Math.min(v.anisotropy, r.getMaxAnisotropy())), n.get(v).__currentAnisotropy = v.anisotropy);
    }
  }
  function Te(E, v) {
    let O = !1;
    E.__webglInit === void 0 && (E.__webglInit = !0, v.addEventListener("dispose", Y));
    const ie = v.source;
    let te = p.get(ie);
    te === void 0 && (te = {}, p.set(ie, te));
    const se = J(v);
    if (se !== E.__cacheKey) {
      te[se] === void 0 && (te[se] = {
        texture: i.createTexture(),
        usedTimes: 0
      }, o.memory.textures++, O = !0), te[se].usedTimes++;
      const Ee = te[E.__cacheKey];
      Ee !== void 0 && (te[E.__cacheKey].usedTimes--, Ee.usedTimes === 0 && H(v)), E.__cacheKey = se, E.__webglTexture = te[se].texture;
    }
    return O;
  }
  function ke(E, v, O) {
    let ie = i.TEXTURE_2D;
    (v.isDataArrayTexture || v.isCompressedArrayTexture) && (ie = i.TEXTURE_2D_ARRAY), v.isData3DTexture && (ie = i.TEXTURE_3D);
    const te = Te(E, v), se = v.source;
    t.bindTexture(ie, E.__webglTexture, i.TEXTURE0 + O);
    const Ee = n.get(se);
    if (se.version !== Ee.__version || te === !0) {
      t.activeTexture(i.TEXTURE0 + O), i.pixelStorei(i.UNPACK_FLIP_Y_WEBGL, v.flipY), i.pixelStorei(i.UNPACK_PREMULTIPLY_ALPHA_WEBGL, v.premultiplyAlpha), i.pixelStorei(i.UNPACK_ALIGNMENT, v.unpackAlignment), i.pixelStorei(i.UNPACK_COLORSPACE_CONVERSION_WEBGL, i.NONE);
      const ae = C(v) && T(v.image) === !1;
      let z = _(v.image, ae, !1, h);
      z = Ze(v, z);
      const R = T(z) || a, K = s.convert(v.format, v.colorSpace);
      let _e = s.convert(v.type), he = V(v.internalFormat, K, _e, v.colorSpace);
      be(ie, v, R);
      let ge;
      const De = v.mipmaps, Ve = a && v.isVideoTexture !== !0, P = Ee.__version === void 0 || te === !0, fe = M(v, z, R);
      if (v.isDepthTexture)
        he = i.DEPTH_COMPONENT, a ? v.type === rn ? he = i.DEPTH_COMPONENT32F : v.type === nn ? he = i.DEPTH_COMPONENT24 : v.type === xn ? he = i.DEPTH24_STENCIL8 : he = i.DEPTH_COMPONENT16 : v.type === rn && console.error("WebGLRenderer: Floating point depth texture requires WebGL2."), v.format === Mn && he === i.DEPTH_COMPONENT && v.type !== as && v.type !== nn && (console.warn("THREE.WebGLRenderer: Use UnsignedShortType or UnsignedIntType for DepthFormat DepthTexture."), v.type = nn, _e = s.convert(v.type)), v.format === ni && he === i.DEPTH_COMPONENT && (he = i.DEPTH_STENCIL, v.type !== xn && (console.warn("THREE.WebGLRenderer: Use UnsignedInt248Type for DepthStencilFormat DepthTexture."), v.type = xn, _e = s.convert(v.type))), P && (Ve ? t.texStorage2D(i.TEXTURE_2D, 1, he, z.width, z.height) : t.texImage2D(i.TEXTURE_2D, 0, he, z.width, z.height, 0, K, _e, null));
      else if (v.isDataTexture)
        if (De.length > 0 && R) {
          Ve && P && t.texStorage2D(i.TEXTURE_2D, fe, he, De[0].width, De[0].height);
          for (let F = 0, ne = De.length; F < ne; F++)
            ge = De[F], Ve ? t.texSubImage2D(i.TEXTURE_2D, F, 0, 0, ge.width, ge.height, K, _e, ge.data) : t.texImage2D(i.TEXTURE_2D, F, he, ge.width, ge.height, 0, K, _e, ge.data);
          v.generateMipmaps = !1;
        } else
          Ve ? (P && t.texStorage2D(i.TEXTURE_2D, fe, he, z.width, z.height), t.texSubImage2D(i.TEXTURE_2D, 0, 0, 0, z.width, z.height, K, _e, z.data)) : t.texImage2D(i.TEXTURE_2D, 0, he, z.width, z.height, 0, K, _e, z.data);
      else if (v.isCompressedTexture)
        if (v.isCompressedArrayTexture) {
          Ve && P && t.texStorage3D(i.TEXTURE_2D_ARRAY, fe, he, De[0].width, De[0].height, z.depth);
          for (let F = 0, ne = De.length; F < ne; F++)
            ge = De[F], v.format !== Ut ? K !== null ? Ve ? t.compressedTexSubImage3D(i.TEXTURE_2D_ARRAY, F, 0, 0, 0, ge.width, ge.height, z.depth, K, ge.data, 0, 0) : t.compressedTexImage3D(i.TEXTURE_2D_ARRAY, F, he, ge.width, ge.height, z.depth, 0, ge.data, 0, 0) : console.warn("THREE.WebGLRenderer: Attempt to load unsupported compressed texture format in .uploadTexture()") : Ve ? t.texSubImage3D(i.TEXTURE_2D_ARRAY, F, 0, 0, 0, ge.width, ge.height, z.depth, K, _e, ge.data) : t.texImage3D(i.TEXTURE_2D_ARRAY, F, he, ge.width, ge.height, z.depth, 0, K, _e, ge.data);
        } else {
          Ve && P && t.texStorage2D(i.TEXTURE_2D, fe, he, De[0].width, De[0].height);
          for (let F = 0, ne = De.length; F < ne; F++)
            ge = De[F], v.format !== Ut ? K !== null ? Ve ? t.compressedTexSubImage2D(i.TEXTURE_2D, F, 0, 0, ge.width, ge.height, K, ge.data) : t.compressedTexImage2D(i.TEXTURE_2D, F, he, ge.width, ge.height, 0, ge.data) : console.warn("THREE.WebGLRenderer: Attempt to load unsupported compressed texture format in .uploadTexture()") : Ve ? t.texSubImage2D(i.TEXTURE_2D, F, 0, 0, ge.width, ge.height, K, _e, ge.data) : t.texImage2D(i.TEXTURE_2D, F, he, ge.width, ge.height, 0, K, _e, ge.data);
        }
      else if (v.isDataArrayTexture)
        Ve ? (P && t.texStorage3D(i.TEXTURE_2D_ARRAY, fe, he, z.width, z.height, z.depth), t.texSubImage3D(i.TEXTURE_2D_ARRAY, 0, 0, 0, 0, z.width, z.height, z.depth, K, _e, z.data)) : t.texImage3D(i.TEXTURE_2D_ARRAY, 0, he, z.width, z.height, z.depth, 0, K, _e, z.data);
      else if (v.isData3DTexture)
        Ve ? (P && t.texStorage3D(i.TEXTURE_3D, fe, he, z.width, z.height, z.depth), t.texSubImage3D(i.TEXTURE_3D, 0, 0, 0, 0, z.width, z.height, z.depth, K, _e, z.data)) : t.texImage3D(i.TEXTURE_3D, 0, he, z.width, z.height, z.depth, 0, K, _e, z.data);
      else if (v.isFramebufferTexture) {
        if (P)
          if (Ve)
            t.texStorage2D(i.TEXTURE_2D, fe, he, z.width, z.height);
          else {
            let F = z.width, ne = z.height;
            for (let de = 0; de < fe; de++)
              t.texImage2D(i.TEXTURE_2D, de, he, F, ne, 0, K, _e, null), F >>= 1, ne >>= 1;
          }
      } else if (De.length > 0 && R) {
        Ve && P && t.texStorage2D(i.TEXTURE_2D, fe, he, De[0].width, De[0].height);
        for (let F = 0, ne = De.length; F < ne; F++)
          ge = De[F], Ve ? t.texSubImage2D(i.TEXTURE_2D, F, 0, 0, K, _e, ge) : t.texImage2D(i.TEXTURE_2D, F, he, K, _e, ge);
        v.generateMipmaps = !1;
      } else
        Ve ? (P && t.texStorage2D(i.TEXTURE_2D, fe, he, z.width, z.height), t.texSubImage2D(i.TEXTURE_2D, 0, 0, 0, K, _e, z)) : t.texImage2D(i.TEXTURE_2D, 0, he, K, _e, z);
      L(v, R) && w(ie), Ee.__version = se.version, v.onUpdate && v.onUpdate(v);
    }
    E.__version = v.version;
  }
  function Ye(E, v, O) {
    if (v.image.length !== 6)
      return;
    const ie = Te(E, v), te = v.source;
    t.bindTexture(i.TEXTURE_CUBE_MAP, E.__webglTexture, i.TEXTURE0 + O);
    const se = n.get(te);
    if (te.version !== se.__version || ie === !0) {
      t.activeTexture(i.TEXTURE0 + O), i.pixelStorei(i.UNPACK_FLIP_Y_WEBGL, v.flipY), i.pixelStorei(i.UNPACK_PREMULTIPLY_ALPHA_WEBGL, v.premultiplyAlpha), i.pixelStorei(i.UNPACK_ALIGNMENT, v.unpackAlignment), i.pixelStorei(i.UNPACK_COLORSPACE_CONVERSION_WEBGL, i.NONE);
      const Ee = v.isCompressedTexture || v.image[0].isCompressedTexture, ae = v.image[0] && v.image[0].isDataTexture, z = [];
      for (let F = 0; F < 6; F++)
        !Ee && !ae ? z[F] = _(v.image[F], !1, !0, c) : z[F] = ae ? v.image[F].image : v.image[F], z[F] = Ze(v, z[F]);
      const R = z[0], K = T(R) || a, _e = s.convert(v.format, v.colorSpace), he = s.convert(v.type), ge = V(v.internalFormat, _e, he, v.colorSpace), De = a && v.isVideoTexture !== !0, Ve = se.__version === void 0 || ie === !0;
      let P = M(v, R, K);
      be(i.TEXTURE_CUBE_MAP, v, K);
      let fe;
      if (Ee) {
        De && Ve && t.texStorage2D(i.TEXTURE_CUBE_MAP, P, ge, R.width, R.height);
        for (let F = 0; F < 6; F++) {
          fe = z[F].mipmaps;
          for (let ne = 0; ne < fe.length; ne++) {
            const de = fe[ne];
            v.format !== Ut ? _e !== null ? De ? t.compressedTexSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne, 0, 0, de.width, de.height, _e, de.data) : t.compressedTexImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne, ge, de.width, de.height, 0, de.data) : console.warn("THREE.WebGLRenderer: Attempt to load unsupported compressed texture format in .setTextureCube()") : De ? t.texSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne, 0, 0, de.width, de.height, _e, he, de.data) : t.texImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne, ge, de.width, de.height, 0, _e, he, de.data);
          }
        }
      } else {
        fe = v.mipmaps, De && Ve && (fe.length > 0 && P++, t.texStorage2D(i.TEXTURE_CUBE_MAP, P, ge, z[0].width, z[0].height));
        for (let F = 0; F < 6; F++)
          if (ae) {
            De ? t.texSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, 0, 0, 0, z[F].width, z[F].height, _e, he, z[F].data) : t.texImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, 0, ge, z[F].width, z[F].height, 0, _e, he, z[F].data);
            for (let ne = 0; ne < fe.length; ne++) {
              const Fe = fe[ne].image[F].image;
              De ? t.texSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne + 1, 0, 0, Fe.width, Fe.height, _e, he, Fe.data) : t.texImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne + 1, ge, Fe.width, Fe.height, 0, _e, he, Fe.data);
            }
          } else {
            De ? t.texSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, 0, 0, 0, _e, he, z[F]) : t.texImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, 0, ge, _e, he, z[F]);
            for (let ne = 0; ne < fe.length; ne++) {
              const de = fe[ne];
              De ? t.texSubImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne + 1, 0, 0, _e, he, de.image[F]) : t.texImage2D(i.TEXTURE_CUBE_MAP_POSITIVE_X + F, ne + 1, ge, _e, he, de.image[F]);
            }
          }
      }
      L(v, K) && w(i.TEXTURE_CUBE_MAP), se.__version = te.version, v.onUpdate && v.onUpdate(v);
    }
    E.__version = v.version;
  }
  function we(E, v, O, ie, te, se) {
    const Ee = s.convert(O.format, O.colorSpace), ae = s.convert(O.type), z = V(O.internalFormat, Ee, ae, O.colorSpace);
    if (!n.get(v).__hasExternalTextures) {
      const K = Math.max(1, v.width >> se), _e = Math.max(1, v.height >> se);
      te === i.TEXTURE_3D || te === i.TEXTURE_2D_ARRAY ? t.texImage3D(te, se, z, K, _e, v.depth, 0, Ee, ae, null) : t.texImage2D(te, se, z, K, _e, 0, Ee, ae, null);
    }
    t.bindFramebuffer(i.FRAMEBUFFER, E), Ce(v) ? u.framebufferTexture2DMultisampleEXT(i.FRAMEBUFFER, ie, te, n.get(O).__webglTexture, 0, Me(v)) : (te === i.TEXTURE_2D || te >= i.TEXTURE_CUBE_MAP_POSITIVE_X && te <= i.TEXTURE_CUBE_MAP_NEGATIVE_Z) && i.framebufferTexture2D(i.FRAMEBUFFER, ie, te, n.get(O).__webglTexture, se), t.bindFramebuffer(i.FRAMEBUFFER, null);
  }
  function b(E, v, O) {
    if (i.bindRenderbuffer(i.RENDERBUFFER, E), v.depthBuffer && !v.stencilBuffer) {
      let ie = i.DEPTH_COMPONENT16;
      if (O || Ce(v)) {
        const te = v.depthTexture;
        te && te.isDepthTexture && (te.type === rn ? ie = i.DEPTH_COMPONENT32F : te.type === nn && (ie = i.DEPTH_COMPONENT24));
        const se = Me(v);
        Ce(v) ? u.renderbufferStorageMultisampleEXT(i.RENDERBUFFER, se, ie, v.width, v.height) : i.renderbufferStorageMultisample(i.RENDERBUFFER, se, ie, v.width, v.height);
      } else
        i.renderbufferStorage(i.RENDERBUFFER, ie, v.width, v.height);
      i.framebufferRenderbuffer(i.FRAMEBUFFER, i.DEPTH_ATTACHMENT, i.RENDERBUFFER, E);
    } else if (v.depthBuffer && v.stencilBuffer) {
      const ie = Me(v);
      O && Ce(v) === !1 ? i.renderbufferStorageMultisample(i.RENDERBUFFER, ie, i.DEPTH24_STENCIL8, v.width, v.height) : Ce(v) ? u.renderbufferStorageMultisampleEXT(i.RENDERBUFFER, ie, i.DEPTH24_STENCIL8, v.width, v.height) : i.renderbufferStorage(i.RENDERBUFFER, i.DEPTH_STENCIL, v.width, v.height), i.framebufferRenderbuffer(i.FRAMEBUFFER, i.DEPTH_STENCIL_ATTACHMENT, i.RENDERBUFFER, E);
    } else {
      const ie = v.isWebGLMultipleRenderTargets === !0 ? v.texture : [v.texture];
      for (let te = 0; te < ie.length; te++) {
        const se = ie[te], Ee = s.convert(se.format, se.colorSpace), ae = s.convert(se.type), z = V(se.internalFormat, Ee, ae, se.colorSpace), R = Me(v);
        O && Ce(v) === !1 ? i.renderbufferStorageMultisample(i.RENDERBUFFER, R, z, v.width, v.height) : Ce(v) ? u.renderbufferStorageMultisampleEXT(i.RENDERBUFFER, R, z, v.width, v.height) : i.renderbufferStorage(i.RENDERBUFFER, z, v.width, v.height);
      }
    }
    i.bindRenderbuffer(i.RENDERBUFFER, null);
  }
  function le(E, v) {
    if (v && v.isWebGLCubeRenderTarget)
      throw new Error("Depth Texture with cube render targets is not supported");
    if (t.bindFramebuffer(i.FRAMEBUFFER, E), !(v.depthTexture && v.depthTexture.isDepthTexture))
      throw new Error("renderTarget.depthTexture must be an instance of THREE.DepthTexture");
    (!n.get(v.depthTexture).__webglTexture || v.depthTexture.image.width !== v.width || v.depthTexture.image.height !== v.height) && (v.depthTexture.image.width = v.width, v.depthTexture.image.height = v.height, v.depthTexture.needsUpdate = !0), ee(v.depthTexture, 0);
    const ie = n.get(v.depthTexture).__webglTexture, te = Me(v);
    if (v.depthTexture.format === Mn)
      Ce(v) ? u.framebufferTexture2DMultisampleEXT(i.FRAMEBUFFER, i.DEPTH_ATTACHMENT, i.TEXTURE_2D, ie, 0, te) : i.framebufferTexture2D(i.FRAMEBUFFER, i.DEPTH_ATTACHMENT, i.TEXTURE_2D, ie, 0);
    else if (v.depthTexture.format === ni)
      Ce(v) ? u.framebufferTexture2DMultisampleEXT(i.FRAMEBUFFER, i.DEPTH_STENCIL_ATTACHMENT, i.TEXTURE_2D, ie, 0, te) : i.framebufferTexture2D(i.FRAMEBUFFER, i.DEPTH_STENCIL_ATTACHMENT, i.TEXTURE_2D, ie, 0);
    else
      throw new Error("Unknown depthTexture format");
  }
  function Z(E) {
    const v = n.get(E), O = E.isWebGLCubeRenderTarget === !0;
    if (E.depthTexture && !v.__autoAllocateDepthBuffer) {
      if (O)
        throw new Error("target.depthTexture not supported in Cube render targets");
      le(v.__webglFramebuffer, E);
    } else if (O) {
      v.__webglDepthbuffer = [];
      for (let ie = 0; ie < 6; ie++)
        t.bindFramebuffer(i.FRAMEBUFFER, v.__webglFramebuffer[ie]), v.__webglDepthbuffer[ie] = i.createRenderbuffer(), b(v.__webglDepthbuffer[ie], E, !1);
    } else
      t.bindFramebuffer(i.FRAMEBUFFER, v.__webglFramebuffer), v.__webglDepthbuffer = i.createRenderbuffer(), b(v.__webglDepthbuffer, E, !1);
    t.bindFramebuffer(i.FRAMEBUFFER, null);
  }
  function re(E, v, O) {
    const ie = n.get(E);
    v !== void 0 && we(ie.__webglFramebuffer, E, E.texture, i.COLOR_ATTACHMENT0, i.TEXTURE_2D, 0), O !== void 0 && Z(E);
  }
  function $(E) {
    const v = E.texture, O = n.get(E), ie = n.get(v);
    E.addEventListener("dispose", ce), E.isWebGLMultipleRenderTargets !== !0 && (ie.__webglTexture === void 0 && (ie.__webglTexture = i.createTexture()), ie.__version = v.version, o.memory.textures++);
    const te = E.isWebGLCubeRenderTarget === !0, se = E.isWebGLMultipleRenderTargets === !0, Ee = T(E) || a;
    if (te) {
      O.__webglFramebuffer = [];
      for (let ae = 0; ae < 6; ae++)
        if (a && v.mipmaps && v.mipmaps.length > 0) {
          O.__webglFramebuffer[ae] = [];
          for (let z = 0; z < v.mipmaps.length; z++)
            O.__webglFramebuffer[ae][z] = i.createFramebuffer();
        } else
          O.__webglFramebuffer[ae] = i.createFramebuffer();
    } else {
      if (a && v.mipmaps && v.mipmaps.length > 0) {
        O.__webglFramebuffer = [];
        for (let ae = 0; ae < v.mipmaps.length; ae++)
          O.__webglFramebuffer[ae] = i.createFramebuffer();
      } else
        O.__webglFramebuffer = i.createFramebuffer();
      if (se)
        if (r.drawBuffers) {
          const ae = E.texture;
          for (let z = 0, R = ae.length; z < R; z++) {
            const K = n.get(ae[z]);
            K.__webglTexture === void 0 && (K.__webglTexture = i.createTexture(), o.memory.textures++);
          }
        } else
          console.warn("THREE.WebGLRenderer: WebGLMultipleRenderTargets can only be used with WebGL2 or WEBGL_draw_buffers extension.");
      if (a && E.samples > 0 && Ce(E) === !1) {
        const ae = se ? v : [v];
        O.__webglMultisampledFramebuffer = i.createFramebuffer(), O.__webglColorRenderbuffer = [], t.bindFramebuffer(i.FRAMEBUFFER, O.__webglMultisampledFramebuffer);
        for (let z = 0; z < ae.length; z++) {
          const R = ae[z];
          O.__webglColorRenderbuffer[z] = i.createRenderbuffer(), i.bindRenderbuffer(i.RENDERBUFFER, O.__webglColorRenderbuffer[z]);
          const K = s.convert(R.format, R.colorSpace), _e = s.convert(R.type), he = V(R.internalFormat, K, _e, R.colorSpace, E.isXRRenderTarget === !0), ge = Me(E);
          i.renderbufferStorageMultisample(i.RENDERBUFFER, ge, he, E.width, E.height), i.framebufferRenderbuffer(i.FRAMEBUFFER, i.COLOR_ATTACHMENT0 + z, i.RENDERBUFFER, O.__webglColorRenderbuffer[z]);
        }
        i.bindRenderbuffer(i.RENDERBUFFER, null), E.depthBuffer && (O.__webglDepthRenderbuffer = i.createRenderbuffer(), b(O.__webglDepthRenderbuffer, E, !0)), t.bindFramebuffer(i.FRAMEBUFFER, null);
      }
    }
    if (te) {
      t.bindTexture(i.TEXTURE_CUBE_MAP, ie.__webglTexture), be(i.TEXTURE_CUBE_MAP, v, Ee);
      for (let ae = 0; ae < 6; ae++)
        if (a && v.mipmaps && v.mipmaps.length > 0)
          for (let z = 0; z < v.mipmaps.length; z++)
            we(O.__webglFramebuffer[ae][z], E, v, i.COLOR_ATTACHMENT0, i.TEXTURE_CUBE_MAP_POSITIVE_X + ae, z);
        else
          we(O.__webglFramebuffer[ae], E, v, i.COLOR_ATTACHMENT0, i.TEXTURE_CUBE_MAP_POSITIVE_X + ae, 0);
      L(v, Ee) && w(i.TEXTURE_CUBE_MAP), t.unbindTexture();
    } else if (se) {
      const ae = E.texture;
      for (let z = 0, R = ae.length; z < R; z++) {
        const K = ae[z], _e = n.get(K);
        t.bindTexture(i.TEXTURE_2D, _e.__webglTexture), be(i.TEXTURE_2D, K, Ee), we(O.__webglFramebuffer, E, K, i.COLOR_ATTACHMENT0 + z, i.TEXTURE_2D, 0), L(K, Ee) && w(i.TEXTURE_2D);
      }
      t.unbindTexture();
    } else {
      let ae = i.TEXTURE_2D;
      if ((E.isWebGL3DRenderTarget || E.isWebGLArrayRenderTarget) && (a ? ae = E.isWebGL3DRenderTarget ? i.TEXTURE_3D : i.TEXTURE_2D_ARRAY : console.error("THREE.WebGLTextures: THREE.Data3DTexture and THREE.DataArrayTexture only supported with WebGL2.")), t.bindTexture(ae, ie.__webglTexture), be(ae, v, Ee), a && v.mipmaps && v.mipmaps.length > 0)
        for (let z = 0; z < v.mipmaps.length; z++)
          we(O.__webglFramebuffer[z], E, v, i.COLOR_ATTACHMENT0, ae, z);
      else
        we(O.__webglFramebuffer, E, v, i.COLOR_ATTACHMENT0, ae, 0);
      L(v, Ee) && w(ae), t.unbindTexture();
    }
    E.depthBuffer && Z(E);
  }
  function ye(E) {
    const v = T(E) || a, O = E.isWebGLMultipleRenderTargets === !0 ? E.texture : [E.texture];
    for (let ie = 0, te = O.length; ie < te; ie++) {
      const se = O[ie];
      if (L(se, v)) {
        const Ee = E.isWebGLCubeRenderTarget ? i.TEXTURE_CUBE_MAP : i.TEXTURE_2D, ae = n.get(se).__webglTexture;
        t.bindTexture(Ee, ae), w(Ee), t.unbindTexture();
      }
    }
  }
  function xe(E) {
    if (a && E.samples > 0 && Ce(E) === !1) {
      const v = E.isWebGLMultipleRenderTargets ? E.texture : [E.texture], O = E.width, ie = E.height;
      let te = i.COLOR_BUFFER_BIT;
      const se = [], Ee = E.stencilBuffer ? i.DEPTH_STENCIL_ATTACHMENT : i.DEPTH_ATTACHMENT, ae = n.get(E), z = E.isWebGLMultipleRenderTargets === !0;
      if (z)
        for (let R = 0; R < v.length; R++)
          t.bindFramebuffer(i.FRAMEBUFFER, ae.__webglMultisampledFramebuffer), i.framebufferRenderbuffer(i.FRAMEBUFFER, i.COLOR_ATTACHMENT0 + R, i.RENDERBUFFER, null), t.bindFramebuffer(i.FRAMEBUFFER, ae.__webglFramebuffer), i.framebufferTexture2D(i.DRAW_FRAMEBUFFER, i.COLOR_ATTACHMENT0 + R, i.TEXTURE_2D, null, 0);
      t.bindFramebuffer(i.READ_FRAMEBUFFER, ae.__webglMultisampledFramebuffer), t.bindFramebuffer(i.DRAW_FRAMEBUFFER, ae.__webglFramebuffer);
      for (let R = 0; R < v.length; R++) {
        se.push(i.COLOR_ATTACHMENT0 + R), E.depthBuffer && se.push(Ee);
        const K = ae.__ignoreDepthValues !== void 0 ? ae.__ignoreDepthValues : !1;
        if (K === !1 && (E.depthBuffer && (te |= i.DEPTH_BUFFER_BIT), E.stencilBuffer && (te |= i.STENCIL_BUFFER_BIT)), z && i.framebufferRenderbuffer(i.READ_FRAMEBUFFER, i.COLOR_ATTACHMENT0, i.RENDERBUFFER, ae.__webglColorRenderbuffer[R]), K === !0 && (i.invalidateFramebuffer(i.READ_FRAMEBUFFER, [Ee]), i.invalidateFramebuffer(i.DRAW_FRAMEBUFFER, [Ee])), z) {
          const _e = n.get(v[R]).__webglTexture;
          i.framebufferTexture2D(i.DRAW_FRAMEBUFFER, i.COLOR_ATTACHMENT0, i.TEXTURE_2D, _e, 0);
        }
        i.blitFramebuffer(0, 0, O, ie, 0, 0, O, ie, te, i.NEAREST), m && i.invalidateFramebuffer(i.READ_FRAMEBUFFER, se);
      }
      if (t.bindFramebuffer(i.READ_FRAMEBUFFER, null), t.bindFramebuffer(i.DRAW_FRAMEBUFFER, null), z)
        for (let R = 0; R < v.length; R++) {
          t.bindFramebuffer(i.FRAMEBUFFER, ae.__webglMultisampledFramebuffer), i.framebufferRenderbuffer(i.FRAMEBUFFER, i.COLOR_ATTACHMENT0 + R, i.RENDERBUFFER, ae.__webglColorRenderbuffer[R]);
          const K = n.get(v[R]).__webglTexture;
          t.bindFramebuffer(i.FRAMEBUFFER, ae.__webglFramebuffer), i.framebufferTexture2D(i.DRAW_FRAMEBUFFER, i.COLOR_ATTACHMENT0 + R, i.TEXTURE_2D, K, 0);
        }
      t.bindFramebuffer(i.DRAW_FRAMEBUFFER, ae.__webglMultisampledFramebuffer);
    }
  }
  function Me(E) {
    return Math.min(f, E.samples);
  }
  function Ce(E) {
    const v = n.get(E);
    return a && E.samples > 0 && e.has("WEBGL_multisampled_render_to_texture") === !0 && v.__useRenderToTexture !== !1;
  }
  function ze(E) {
    const v = o.render.frame;
    g.get(E) !== v && (g.set(E, v), E.update());
  }
  function Ze(E, v) {
    const O = E.colorSpace, ie = E.format, te = E.type;
    return E.isCompressedTexture === !0 || E.format === Jr || O !== Ft && O !== En && (O === Oe ? a === !1 ? e.has("EXT_sRGB") === !0 && ie === Ut ? (E.format = Jr, E.minFilter = Tt, E.generateMipmaps = !1) : v = po.sRGBToLinear(v) : (ie !== Ut || te !== on) && console.warn("THREE.WebGLTextures: sRGB encoded textures have to use RGBAFormat and UnsignedByteType.") : console.error("THREE.WebGLTextures: Unsupported texture color space:", O)), v;
  }
  this.allocateTextureUnit = j, this.resetTextureUnits = X, this.setTexture2D = ee, this.setTexture2DArray = I, this.setTexture3D = q, this.setTextureCube = pe, this.rebindTextures = re, this.setupRenderTarget = $, this.updateRenderTargetMipmap = ye, this.updateMultisampleRenderTarget = xe, this.setupDepthRenderbuffer = Z, this.setupFrameBufferTexture = we, this.useMultisampledRTT = Ce;
}
function ap(i, e, t) {
  const n = t.isWebGL2;
  function r(s, o = En) {
    let a;
    if (s === on)
      return i.UNSIGNED_BYTE;
    if (s === ro)
      return i.UNSIGNED_SHORT_4_4_4_4;
    if (s === so)
      return i.UNSIGNED_SHORT_5_5_5_1;
    if (s === Al)
      return i.BYTE;
    if (s === bl)
      return i.SHORT;
    if (s === as)
      return i.UNSIGNED_SHORT;
    if (s === io)
      return i.INT;
    if (s === nn)
      return i.UNSIGNED_INT;
    if (s === rn)
      return i.FLOAT;
    if (s === xi)
      return n ? i.HALF_FLOAT : (a = e.get("OES_texture_half_float"), a !== null ? a.HALF_FLOAT_OES : null);
    if (s === wl)
      return i.ALPHA;
    if (s === Ut)
      return i.RGBA;
    if (s === Rl)
      return i.LUMINANCE;
    if (s === Cl)
      return i.LUMINANCE_ALPHA;
    if (s === Mn)
      return i.DEPTH_COMPONENT;
    if (s === ni)
      return i.DEPTH_STENCIL;
    if (s === Jr)
      return a = e.get("EXT_sRGB"), a !== null ? a.SRGB_ALPHA_EXT : null;
    if (s === Pl)
      return i.RED;
    if (s === ao)
      return i.RED_INTEGER;
    if (s === Ll)
      return i.RG;
    if (s === oo)
      return i.RG_INTEGER;
    if (s === lo)
      return i.RGBA_INTEGER;
    if (s === pr || s === mr || s === gr || s === _r)
      if (o === Oe)
        if (a = e.get("WEBGL_compressed_texture_s3tc_srgb"), a !== null) {
          if (s === pr)
            return a.COMPRESSED_SRGB_S3TC_DXT1_EXT;
          if (s === mr)
            return a.COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT;
          if (s === gr)
            return a.COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT;
          if (s === _r)
            return a.COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT;
        } else
          return null;
      else if (a = e.get("WEBGL_compressed_texture_s3tc"), a !== null) {
        if (s === pr)
          return a.COMPRESSED_RGB_S3TC_DXT1_EXT;
        if (s === mr)
          return a.COMPRESSED_RGBA_S3TC_DXT1_EXT;
        if (s === gr)
          return a.COMPRESSED_RGBA_S3TC_DXT3_EXT;
        if (s === _r)
          return a.COMPRESSED_RGBA_S3TC_DXT5_EXT;
      } else
        return null;
    if (s === Ps || s === Ls || s === Us || s === Ds)
      if (a = e.get("WEBGL_compressed_texture_pvrtc"), a !== null) {
        if (s === Ps)
          return a.COMPRESSED_RGB_PVRTC_4BPPV1_IMG;
        if (s === Ls)
          return a.COMPRESSED_RGB_PVRTC_2BPPV1_IMG;
        if (s === Us)
          return a.COMPRESSED_RGBA_PVRTC_4BPPV1_IMG;
        if (s === Ds)
          return a.COMPRESSED_RGBA_PVRTC_2BPPV1_IMG;
      } else
        return null;
    if (s === Ul)
      return a = e.get("WEBGL_compressed_texture_etc1"), a !== null ? a.COMPRESSED_RGB_ETC1_WEBGL : null;
    if (s === Is || s === Ns)
      if (a = e.get("WEBGL_compressed_texture_etc"), a !== null) {
        if (s === Is)
          return o === Oe ? a.COMPRESSED_SRGB8_ETC2 : a.COMPRESSED_RGB8_ETC2;
        if (s === Ns)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ETC2_EAC : a.COMPRESSED_RGBA8_ETC2_EAC;
      } else
        return null;
    if (s === Os || s === Fs || s === Bs || s === zs || s === Hs || s === Gs || s === Vs || s === ks || s === Ws || s === Xs || s === Ys || s === qs || s === js || s === Zs)
      if (a = e.get("WEBGL_compressed_texture_astc"), a !== null) {
        if (s === Os)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_4x4_KHR : a.COMPRESSED_RGBA_ASTC_4x4_KHR;
        if (s === Fs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_5x4_KHR : a.COMPRESSED_RGBA_ASTC_5x4_KHR;
        if (s === Bs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_5x5_KHR : a.COMPRESSED_RGBA_ASTC_5x5_KHR;
        if (s === zs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_6x5_KHR : a.COMPRESSED_RGBA_ASTC_6x5_KHR;
        if (s === Hs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_6x6_KHR : a.COMPRESSED_RGBA_ASTC_6x6_KHR;
        if (s === Gs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_8x5_KHR : a.COMPRESSED_RGBA_ASTC_8x5_KHR;
        if (s === Vs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_8x6_KHR : a.COMPRESSED_RGBA_ASTC_8x6_KHR;
        if (s === ks)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_8x8_KHR : a.COMPRESSED_RGBA_ASTC_8x8_KHR;
        if (s === Ws)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_10x5_KHR : a.COMPRESSED_RGBA_ASTC_10x5_KHR;
        if (s === Xs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_10x6_KHR : a.COMPRESSED_RGBA_ASTC_10x6_KHR;
        if (s === Ys)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_10x8_KHR : a.COMPRESSED_RGBA_ASTC_10x8_KHR;
        if (s === qs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_10x10_KHR : a.COMPRESSED_RGBA_ASTC_10x10_KHR;
        if (s === js)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_12x10_KHR : a.COMPRESSED_RGBA_ASTC_12x10_KHR;
        if (s === Zs)
          return o === Oe ? a.COMPRESSED_SRGB8_ALPHA8_ASTC_12x12_KHR : a.COMPRESSED_RGBA_ASTC_12x12_KHR;
      } else
        return null;
    if (s === vr)
      if (a = e.get("EXT_texture_compression_bptc"), a !== null) {
        if (s === vr)
          return o === Oe ? a.COMPRESSED_SRGB_ALPHA_BPTC_UNORM_EXT : a.COMPRESSED_RGBA_BPTC_UNORM_EXT;
      } else
        return null;
    if (s === Dl || s === Ks || s === Js || s === $s)
      if (a = e.get("EXT_texture_compression_rgtc"), a !== null) {
        if (s === vr)
          return a.COMPRESSED_RED_RGTC1_EXT;
        if (s === Ks)
          return a.COMPRESSED_SIGNED_RED_RGTC1_EXT;
        if (s === Js)
          return a.COMPRESSED_RED_GREEN_RGTC2_EXT;
        if (s === $s)
          return a.COMPRESSED_SIGNED_RED_GREEN_RGTC2_EXT;
      } else
        return null;
    return s === xn ? n ? i.UNSIGNED_INT_24_8 : (a = e.get("WEBGL_depth_texture"), a !== null ? a.UNSIGNED_INT_24_8_WEBGL : null) : i[s] !== void 0 ? i[s] : null;
  }
  return { convert: r };
}
class op extends At {
  constructor(e = []) {
    super(), this.isArrayCamera = !0, this.cameras = e;
  }
}
class Zi extends ht {
  constructor() {
    super(), this.isGroup = !0, this.type = "Group";
  }
}
const lp = { type: "move" };
class Hr {
  constructor() {
    this._targetRay = null, this._grip = null, this._hand = null;
  }
  getHandSpace() {
    return this._hand === null && (this._hand = new Zi(), this._hand.matrixAutoUpdate = !1, this._hand.visible = !1, this._hand.joints = {}, this._hand.inputState = { pinching: !1 }), this._hand;
  }
  getTargetRaySpace() {
    return this._targetRay === null && (this._targetRay = new Zi(), this._targetRay.matrixAutoUpdate = !1, this._targetRay.visible = !1, this._targetRay.hasLinearVelocity = !1, this._targetRay.linearVelocity = new U(), this._targetRay.hasAngularVelocity = !1, this._targetRay.angularVelocity = new U()), this._targetRay;
  }
  getGripSpace() {
    return this._grip === null && (this._grip = new Zi(), this._grip.matrixAutoUpdate = !1, this._grip.visible = !1, this._grip.hasLinearVelocity = !1, this._grip.linearVelocity = new U(), this._grip.hasAngularVelocity = !1, this._grip.angularVelocity = new U()), this._grip;
  }
  dispatchEvent(e) {
    return this._targetRay !== null && this._targetRay.dispatchEvent(e), this._grip !== null && this._grip.dispatchEvent(e), this._hand !== null && this._hand.dispatchEvent(e), this;
  }
  connect(e) {
    if (e && e.hand) {
      const t = this._hand;
      if (t)
        for (const n of e.hand.values())
          this._getHandJoint(t, n);
    }
    return this.dispatchEvent({ type: "connected", data: e }), this;
  }
  disconnect(e) {
    return this.dispatchEvent({ type: "disconnected", data: e }), this._targetRay !== null && (this._targetRay.visible = !1), this._grip !== null && (this._grip.visible = !1), this._hand !== null && (this._hand.visible = !1), this;
  }
  update(e, t, n) {
    let r = null, s = null, o = null;
    const a = this._targetRay, l = this._grip, c = this._hand;
    if (e && t.session.visibilityState !== "visible-blurred") {
      if (c && e.hand) {
        o = !0;
        for (const x of e.hand.values()) {
          const p = t.getJointPose(x, n), d = this._getHandJoint(c, x);
          p !== null && (d.matrix.fromArray(p.transform.matrix), d.matrix.decompose(d.position, d.rotation, d.scale), d.matrixWorldNeedsUpdate = !0, d.jointRadius = p.radius), d.visible = p !== null;
        }
        const h = c.joints["index-finger-tip"], f = c.joints["thumb-tip"], u = h.position.distanceTo(f.position), m = 0.02, g = 5e-3;
        c.inputState.pinching && u > m + g ? (c.inputState.pinching = !1, this.dispatchEvent({
          type: "pinchend",
          handedness: e.handedness,
          target: this
        })) : !c.inputState.pinching && u <= m - g && (c.inputState.pinching = !0, this.dispatchEvent({
          type: "pinchstart",
          handedness: e.handedness,
          target: this
        }));
      } else
        l !== null && e.gripSpace && (s = t.getPose(e.gripSpace, n), s !== null && (l.matrix.fromArray(s.transform.matrix), l.matrix.decompose(l.position, l.rotation, l.scale), l.matrixWorldNeedsUpdate = !0, s.linearVelocity ? (l.hasLinearVelocity = !0, l.linearVelocity.copy(s.linearVelocity)) : l.hasLinearVelocity = !1, s.angularVelocity ? (l.hasAngularVelocity = !0, l.angularVelocity.copy(s.angularVelocity)) : l.hasAngularVelocity = !1));
      a !== null && (r = t.getPose(e.targetRaySpace, n), r === null && s !== null && (r = s), r !== null && (a.matrix.fromArray(r.transform.matrix), a.matrix.decompose(a.position, a.rotation, a.scale), a.matrixWorldNeedsUpdate = !0, r.linearVelocity ? (a.hasLinearVelocity = !0, a.linearVelocity.copy(r.linearVelocity)) : a.hasLinearVelocity = !1, r.angularVelocity ? (a.hasAngularVelocity = !0, a.angularVelocity.copy(r.angularVelocity)) : a.hasAngularVelocity = !1, this.dispatchEvent(lp)));
    }
    return a !== null && (a.visible = r !== null), l !== null && (l.visible = s !== null), c !== null && (c.visible = o !== null), this;
  }
  // private method
  _getHandJoint(e, t) {
    if (e.joints[t.jointName] === void 0) {
      const n = new Zi();
      n.matrixAutoUpdate = !1, n.visible = !1, e.joints[t.jointName] = n, e.add(n);
    }
    return e.joints[t.jointName];
  }
}
class cp extends St {
  constructor(e, t, n, r, s, o, a, l, c, h) {
    if (h = h !== void 0 ? h : Mn, h !== Mn && h !== ni)
      throw new Error("DepthTexture format must be either THREE.DepthFormat or THREE.DepthStencilFormat");
    n === void 0 && h === Mn && (n = nn), n === void 0 && h === ni && (n = xn), super(null, r, s, o, a, l, h, n, c), this.isDepthTexture = !0, this.image = { width: e, height: t }, this.magFilter = a !== void 0 ? a : mt, this.minFilter = l !== void 0 ? l : mt, this.flipY = !1, this.generateMipmaps = !1, this.compareFunction = null;
  }
  copy(e) {
    return super.copy(e), this.compareFunction = e.compareFunction, this;
  }
  toJSON(e) {
    const t = super.toJSON(e);
    return this.compareFunction !== null && (t.compareFunction = this.compareFunction), t;
  }
}
class hp extends Rn {
  constructor(e, t) {
    super();
    const n = this;
    let r = null, s = 1, o = null, a = "local-floor", l = 1, c = null, h = null, f = null, u = null, m = null, g = null;
    const x = t.getContextAttributes();
    let p = null, d = null;
    const A = [], _ = [], T = new At();
    T.layers.enable(1), T.viewport = new ot();
    const C = new At();
    C.layers.enable(2), C.viewport = new ot();
    const L = [T, C], w = new op();
    w.layers.enable(1), w.layers.enable(2);
    let V = null, M = null;
    this.cameraAutoUpdate = !0, this.enabled = !1, this.isPresenting = !1, this.getController = function(I) {
      let q = A[I];
      return q === void 0 && (q = new Hr(), A[I] = q), q.getTargetRaySpace();
    }, this.getControllerGrip = function(I) {
      let q = A[I];
      return q === void 0 && (q = new Hr(), A[I] = q), q.getGripSpace();
    }, this.getHand = function(I) {
      let q = A[I];
      return q === void 0 && (q = new Hr(), A[I] = q), q.getHandSpace();
    };
    function y(I) {
      const q = _.indexOf(I.inputSource);
      if (q === -1)
        return;
      const pe = A[q];
      pe !== void 0 && (pe.update(I.inputSource, I.frame, c || o), pe.dispatchEvent({ type: I.type, data: I.inputSource }));
    }
    function Y() {
      r.removeEventListener("select", y), r.removeEventListener("selectstart", y), r.removeEventListener("selectend", y), r.removeEventListener("squeeze", y), r.removeEventListener("squeezestart", y), r.removeEventListener("squeezeend", y), r.removeEventListener("end", Y), r.removeEventListener("inputsourceschange", ce);
      for (let I = 0; I < A.length; I++) {
        const q = _[I];
        q !== null && (_[I] = null, A[I].disconnect(q));
      }
      V = null, M = null, e.setRenderTarget(p), m = null, u = null, f = null, r = null, d = null, ee.stop(), n.isPresenting = !1, n.dispatchEvent({ type: "sessionend" });
    }
    this.setFramebufferScaleFactor = function(I) {
      s = I, n.isPresenting === !0 && console.warn("THREE.WebXRManager: Cannot change framebuffer scale while presenting.");
    }, this.setReferenceSpaceType = function(I) {
      a = I, n.isPresenting === !0 && console.warn("THREE.WebXRManager: Cannot change reference space type while presenting.");
    }, this.getReferenceSpace = function() {
      return c || o;
    }, this.setReferenceSpace = function(I) {
      c = I;
    }, this.getBaseLayer = function() {
      return u !== null ? u : m;
    }, this.getBinding = function() {
      return f;
    }, this.getFrame = function() {
      return g;
    }, this.getSession = function() {
      return r;
    }, this.setSession = async function(I) {
      if (r = I, r !== null) {
        if (p = e.getRenderTarget(), r.addEventListener("select", y), r.addEventListener("selectstart", y), r.addEventListener("selectend", y), r.addEventListener("squeeze", y), r.addEventListener("squeezestart", y), r.addEventListener("squeezeend", y), r.addEventListener("end", Y), r.addEventListener("inputsourceschange", ce), x.xrCompatible !== !0 && await t.makeXRCompatible(), r.renderState.layers === void 0 || e.capabilities.isWebGL2 === !1) {
          const q = {
            antialias: r.renderState.layers === void 0 ? x.antialias : !0,
            alpha: !0,
            depth: x.depth,
            stencil: x.stencil,
            framebufferScaleFactor: s
          };
          m = new XRWebGLLayer(r, t, q), r.updateRenderState({ baseLayer: m }), d = new Tn(
            m.framebufferWidth,
            m.framebufferHeight,
            {
              format: Ut,
              type: on,
              colorSpace: e.outputColorSpace,
              stencilBuffer: x.stencil
            }
          );
        } else {
          let q = null, pe = null, me = null;
          x.depth && (me = x.stencil ? t.DEPTH24_STENCIL8 : t.DEPTH_COMPONENT24, q = x.stencil ? ni : Mn, pe = x.stencil ? xn : nn);
          const ve = {
            colorFormat: t.RGBA8,
            depthFormat: me,
            scaleFactor: s
          };
          f = new XRWebGLBinding(r, t), u = f.createProjectionLayer(ve), r.updateRenderState({ layers: [u] }), d = new Tn(
            u.textureWidth,
            u.textureHeight,
            {
              format: Ut,
              type: on,
              depthTexture: new cp(u.textureWidth, u.textureHeight, pe, void 0, void 0, void 0, void 0, void 0, void 0, q),
              stencilBuffer: x.stencil,
              colorSpace: e.outputColorSpace,
              samples: x.antialias ? 4 : 0
            }
          );
          const Ae = e.properties.get(d);
          Ae.__ignoreDepthValues = u.ignoreDepthValues;
        }
        d.isXRRenderTarget = !0, this.setFoveation(l), c = null, o = await r.requestReferenceSpace(a), ee.setContext(r), ee.start(), n.isPresenting = !0, n.dispatchEvent({ type: "sessionstart" });
      }
    }, this.getEnvironmentBlendMode = function() {
      if (r !== null)
        return r.environmentBlendMode;
    };
    function ce(I) {
      for (let q = 0; q < I.removed.length; q++) {
        const pe = I.removed[q], me = _.indexOf(pe);
        me >= 0 && (_[me] = null, A[me].disconnect(pe));
      }
      for (let q = 0; q < I.added.length; q++) {
        const pe = I.added[q];
        let me = _.indexOf(pe);
        if (me === -1) {
          for (let Ae = 0; Ae < A.length; Ae++)
            if (Ae >= _.length) {
              _.push(pe), me = Ae;
              break;
            } else if (_[Ae] === null) {
              _[Ae] = pe, me = Ae;
              break;
            }
          if (me === -1)
            break;
        }
        const ve = A[me];
        ve && ve.connect(pe);
      }
    }
    const B = new U(), H = new U();
    function G(I, q, pe) {
      B.setFromMatrixPosition(q.matrixWorld), H.setFromMatrixPosition(pe.matrixWorld);
      const me = B.distanceTo(H), ve = q.projectionMatrix.elements, Ae = pe.projectionMatrix.elements, be = ve[14] / (ve[10] - 1), Te = ve[14] / (ve[10] + 1), ke = (ve[9] + 1) / ve[5], Ye = (ve[9] - 1) / ve[5], we = (ve[8] - 1) / ve[0], b = (Ae[8] + 1) / Ae[0], le = be * we, Z = be * b, re = me / (-we + b), $ = re * -we;
      q.matrixWorld.decompose(I.position, I.quaternion, I.scale), I.translateX($), I.translateZ(re), I.matrixWorld.compose(I.position, I.quaternion, I.scale), I.matrixWorldInverse.copy(I.matrixWorld).invert();
      const ye = be + re, xe = Te + re, Me = le - $, Ce = Z + (me - $), ze = ke * Te / xe * ye, Ze = Ye * Te / xe * ye;
      I.projectionMatrix.makePerspective(Me, Ce, ze, Ze, ye, xe), I.projectionMatrixInverse.copy(I.projectionMatrix).invert();
    }
    function Q(I, q) {
      q === null ? I.matrixWorld.copy(I.matrix) : I.matrixWorld.multiplyMatrices(q.matrixWorld, I.matrix), I.matrixWorldInverse.copy(I.matrixWorld).invert();
    }
    this.updateCamera = function(I) {
      if (r === null)
        return;
      w.near = C.near = T.near = I.near, w.far = C.far = T.far = I.far, (V !== w.near || M !== w.far) && (r.updateRenderState({
        depthNear: w.near,
        depthFar: w.far
      }), V = w.near, M = w.far);
      const q = I.parent, pe = w.cameras;
      Q(w, q);
      for (let me = 0; me < pe.length; me++)
        Q(pe[me], q);
      pe.length === 2 ? G(w, T, C) : w.projectionMatrix.copy(T.projectionMatrix), X(I, w, q);
    };
    function X(I, q, pe) {
      pe === null ? I.matrix.copy(q.matrixWorld) : (I.matrix.copy(pe.matrixWorld), I.matrix.invert(), I.matrix.multiply(q.matrixWorld)), I.matrix.decompose(I.position, I.quaternion, I.scale), I.updateMatrixWorld(!0);
      const me = I.children;
      for (let ve = 0, Ae = me.length; ve < Ae; ve++)
        me[ve].updateMatrixWorld(!0);
      I.projectionMatrix.copy(q.projectionMatrix), I.projectionMatrixInverse.copy(q.projectionMatrixInverse), I.isPerspectiveCamera && (I.fov = Mi * 2 * Math.atan(1 / I.projectionMatrix.elements[5]), I.zoom = 1);
    }
    this.getCamera = function() {
      return w;
    }, this.getFoveation = function() {
      if (!(u === null && m === null))
        return l;
    }, this.setFoveation = function(I) {
      l = I, u !== null && (u.fixedFoveation = I), m !== null && m.fixedFoveation !== void 0 && (m.fixedFoveation = I);
    };
    let j = null;
    function J(I, q) {
      if (h = q.getViewerPose(c || o), g = q, h !== null) {
        const pe = h.views;
        m !== null && (e.setRenderTargetFramebuffer(d, m.framebuffer), e.setRenderTarget(d));
        let me = !1;
        pe.length !== w.cameras.length && (w.cameras.length = 0, me = !0);
        for (let ve = 0; ve < pe.length; ve++) {
          const Ae = pe[ve];
          let be = null;
          if (m !== null)
            be = m.getViewport(Ae);
          else {
            const ke = f.getViewSubImage(u, Ae);
            be = ke.viewport, ve === 0 && (e.setRenderTargetTextures(
              d,
              ke.colorTexture,
              u.ignoreDepthValues ? void 0 : ke.depthStencilTexture
            ), e.setRenderTarget(d));
          }
          let Te = L[ve];
          Te === void 0 && (Te = new At(), Te.layers.enable(ve), Te.viewport = new ot(), L[ve] = Te), Te.matrix.fromArray(Ae.transform.matrix), Te.matrix.decompose(Te.position, Te.quaternion, Te.scale), Te.projectionMatrix.fromArray(Ae.projectionMatrix), Te.projectionMatrixInverse.copy(Te.projectionMatrix).invert(), Te.viewport.set(be.x, be.y, be.width, be.height), ve === 0 && (w.matrix.copy(Te.matrix), w.matrix.decompose(w.position, w.quaternion, w.scale)), me === !0 && w.cameras.push(Te);
        }
      }
      for (let pe = 0; pe < A.length; pe++) {
        const me = _[pe], ve = A[pe];
        me !== null && ve !== void 0 && ve.update(me, q, c || o);
      }
      j && j(I, q), q.detectedPlanes && n.dispatchEvent({ type: "planesdetected", data: q }), g = null;
    }
    const ee = new To();
    ee.setAnimationLoop(J), this.setAnimationLoop = function(I) {
      j = I;
    }, this.dispose = function() {
    };
  }
}
function up(i, e) {
  function t(p, d) {
    p.matrixAutoUpdate === !0 && p.updateMatrix(), d.value.copy(p.matrix);
  }
  function n(p, d) {
    d.color.getRGB(p.fogColor.value, So(i)), d.isFog ? (p.fogNear.value = d.near, p.fogFar.value = d.far) : d.isFogExp2 && (p.fogDensity.value = d.density);
  }
  function r(p, d, A, _, T) {
    d.isMeshBasicMaterial || d.isMeshLambertMaterial ? s(p, d) : d.isMeshToonMaterial ? (s(p, d), f(p, d)) : d.isMeshPhongMaterial ? (s(p, d), h(p, d)) : d.isMeshStandardMaterial ? (s(p, d), u(p, d), d.isMeshPhysicalMaterial && m(p, d, T)) : d.isMeshMatcapMaterial ? (s(p, d), g(p, d)) : d.isMeshDepthMaterial ? s(p, d) : d.isMeshDistanceMaterial ? (s(p, d), x(p, d)) : d.isMeshNormalMaterial ? s(p, d) : d.isLineBasicMaterial ? (o(p, d), d.isLineDashedMaterial && a(p, d)) : d.isPointsMaterial ? l(p, d, A, _) : d.isSpriteMaterial ? c(p, d) : d.isShadowMaterial ? (p.color.value.copy(d.color), p.opacity.value = d.opacity) : d.isShaderMaterial && (d.uniformsNeedUpdate = !1);
  }
  function s(p, d) {
    p.opacity.value = d.opacity, d.color && p.diffuse.value.copy(d.color), d.emissive && p.emissive.value.copy(d.emissive).multiplyScalar(d.emissiveIntensity), d.map && (p.map.value = d.map, t(d.map, p.mapTransform)), d.alphaMap && (p.alphaMap.value = d.alphaMap, t(d.alphaMap, p.alphaMapTransform)), d.bumpMap && (p.bumpMap.value = d.bumpMap, t(d.bumpMap, p.bumpMapTransform), p.bumpScale.value = d.bumpScale, d.side === gt && (p.bumpScale.value *= -1)), d.normalMap && (p.normalMap.value = d.normalMap, t(d.normalMap, p.normalMapTransform), p.normalScale.value.copy(d.normalScale), d.side === gt && p.normalScale.value.negate()), d.displacementMap && (p.displacementMap.value = d.displacementMap, t(d.displacementMap, p.displacementMapTransform), p.displacementScale.value = d.displacementScale, p.displacementBias.value = d.displacementBias), d.emissiveMap && (p.emissiveMap.value = d.emissiveMap, t(d.emissiveMap, p.emissiveMapTransform)), d.specularMap && (p.specularMap.value = d.specularMap, t(d.specularMap, p.specularMapTransform)), d.alphaTest > 0 && (p.alphaTest.value = d.alphaTest);
    const A = e.get(d).envMap;
    if (A && (p.envMap.value = A, p.flipEnvMap.value = A.isCubeTexture && A.isRenderTargetTexture === !1 ? -1 : 1, p.reflectivity.value = d.reflectivity, p.ior.value = d.ior, p.refractionRatio.value = d.refractionRatio), d.lightMap) {
      p.lightMap.value = d.lightMap;
      const _ = i._useLegacyLights === !0 ? Math.PI : 1;
      p.lightMapIntensity.value = d.lightMapIntensity * _, t(d.lightMap, p.lightMapTransform);
    }
    d.aoMap && (p.aoMap.value = d.aoMap, p.aoMapIntensity.value = d.aoMapIntensity, t(d.aoMap, p.aoMapTransform));
  }
  function o(p, d) {
    p.diffuse.value.copy(d.color), p.opacity.value = d.opacity, d.map && (p.map.value = d.map, t(d.map, p.mapTransform));
  }
  function a(p, d) {
    p.dashSize.value = d.dashSize, p.totalSize.value = d.dashSize + d.gapSize, p.scale.value = d.scale;
  }
  function l(p, d, A, _) {
    p.diffuse.value.copy(d.color), p.opacity.value = d.opacity, p.size.value = d.size * A, p.scale.value = _ * 0.5, d.map && (p.map.value = d.map, t(d.map, p.uvTransform)), d.alphaMap && (p.alphaMap.value = d.alphaMap, t(d.alphaMap, p.alphaMapTransform)), d.alphaTest > 0 && (p.alphaTest.value = d.alphaTest);
  }
  function c(p, d) {
    p.diffuse.value.copy(d.color), p.opacity.value = d.opacity, p.rotation.value = d.rotation, d.map && (p.map.value = d.map, t(d.map, p.mapTransform)), d.alphaMap && (p.alphaMap.value = d.alphaMap, t(d.alphaMap, p.alphaMapTransform)), d.alphaTest > 0 && (p.alphaTest.value = d.alphaTest);
  }
  function h(p, d) {
    p.specular.value.copy(d.specular), p.shininess.value = Math.max(d.shininess, 1e-4);
  }
  function f(p, d) {
    d.gradientMap && (p.gradientMap.value = d.gradientMap);
  }
  function u(p, d) {
    p.metalness.value = d.metalness, d.metalnessMap && (p.metalnessMap.value = d.metalnessMap, t(d.metalnessMap, p.metalnessMapTransform)), p.roughness.value = d.roughness, d.roughnessMap && (p.roughnessMap.value = d.roughnessMap, t(d.roughnessMap, p.roughnessMapTransform)), e.get(d).envMap && (p.envMapIntensity.value = d.envMapIntensity);
  }
  function m(p, d, A) {
    p.ior.value = d.ior, d.sheen > 0 && (p.sheenColor.value.copy(d.sheenColor).multiplyScalar(d.sheen), p.sheenRoughness.value = d.sheenRoughness, d.sheenColorMap && (p.sheenColorMap.value = d.sheenColorMap, t(d.sheenColorMap, p.sheenColorMapTransform)), d.sheenRoughnessMap && (p.sheenRoughnessMap.value = d.sheenRoughnessMap, t(d.sheenRoughnessMap, p.sheenRoughnessMapTransform))), d.clearcoat > 0 && (p.clearcoat.value = d.clearcoat, p.clearcoatRoughness.value = d.clearcoatRoughness, d.clearcoatMap && (p.clearcoatMap.value = d.clearcoatMap, t(d.clearcoatMap, p.clearcoatMapTransform)), d.clearcoatRoughnessMap && (p.clearcoatRoughnessMap.value = d.clearcoatRoughnessMap, t(d.clearcoatRoughnessMap, p.clearcoatRoughnessMapTransform)), d.clearcoatNormalMap && (p.clearcoatNormalMap.value = d.clearcoatNormalMap, t(d.clearcoatNormalMap, p.clearcoatNormalMapTransform), p.clearcoatNormalScale.value.copy(d.clearcoatNormalScale), d.side === gt && p.clearcoatNormalScale.value.negate())), d.iridescence > 0 && (p.iridescence.value = d.iridescence, p.iridescenceIOR.value = d.iridescenceIOR, p.iridescenceThicknessMinimum.value = d.iridescenceThicknessRange[0], p.iridescenceThicknessMaximum.value = d.iridescenceThicknessRange[1], d.iridescenceMap && (p.iridescenceMap.value = d.iridescenceMap, t(d.iridescenceMap, p.iridescenceMapTransform)), d.iridescenceThicknessMap && (p.iridescenceThicknessMap.value = d.iridescenceThicknessMap, t(d.iridescenceThicknessMap, p.iridescenceThicknessMapTransform))), d.transmission > 0 && (p.transmission.value = d.transmission, p.transmissionSamplerMap.value = A.texture, p.transmissionSamplerSize.value.set(A.width, A.height), d.transmissionMap && (p.transmissionMap.value = d.transmissionMap, t(d.transmissionMap, p.transmissionMapTransform)), p.thickness.value = d.thickness, d.thicknessMap && (p.thicknessMap.value = d.thicknessMap, t(d.thicknessMap, p.thicknessMapTransform)), p.attenuationDistance.value = d.attenuationDistance, p.attenuationColor.value.copy(d.attenuationColor)), d.anisotropy > 0 && (p.anisotropyVector.value.set(d.anisotropy * Math.cos(d.anisotropyRotation), d.anisotropy * Math.sin(d.anisotropyRotation)), d.anisotropyMap && (p.anisotropyMap.value = d.anisotropyMap, t(d.anisotropyMap, p.anisotropyMapTransform))), p.specularIntensity.value = d.specularIntensity, p.specularColor.value.copy(d.specularColor), d.specularColorMap && (p.specularColorMap.value = d.specularColorMap, t(d.specularColorMap, p.specularColorMapTransform)), d.specularIntensityMap && (p.specularIntensityMap.value = d.specularIntensityMap, t(d.specularIntensityMap, p.specularIntensityMapTransform));
  }
  function g(p, d) {
    d.matcap && (p.matcap.value = d.matcap);
  }
  function x(p, d) {
    const A = e.get(d).light;
    p.referencePosition.value.setFromMatrixPosition(A.matrixWorld), p.nearDistance.value = A.shadow.camera.near, p.farDistance.value = A.shadow.camera.far;
  }
  return {
    refreshFogUniforms: n,
    refreshMaterialUniforms: r
  };
}
function fp(i, e, t, n) {
  let r = {}, s = {}, o = [];
  const a = t.isWebGL2 ? i.getParameter(i.MAX_UNIFORM_BUFFER_BINDINGS) : 0;
  function l(A, _) {
    const T = _.program;
    n.uniformBlockBinding(A, T);
  }
  function c(A, _) {
    let T = r[A.id];
    T === void 0 && (g(A), T = h(A), r[A.id] = T, A.addEventListener("dispose", p));
    const C = _.program;
    n.updateUBOMapping(A, C);
    const L = e.render.frame;
    s[A.id] !== L && (u(A), s[A.id] = L);
  }
  function h(A) {
    const _ = f();
    A.__bindingPointIndex = _;
    const T = i.createBuffer(), C = A.__size, L = A.usage;
    return i.bindBuffer(i.UNIFORM_BUFFER, T), i.bufferData(i.UNIFORM_BUFFER, C, L), i.bindBuffer(i.UNIFORM_BUFFER, null), i.bindBufferBase(i.UNIFORM_BUFFER, _, T), T;
  }
  function f() {
    for (let A = 0; A < a; A++)
      if (o.indexOf(A) === -1)
        return o.push(A), A;
    return console.error("THREE.WebGLRenderer: Maximum number of simultaneously usable uniforms groups reached."), 0;
  }
  function u(A) {
    const _ = r[A.id], T = A.uniforms, C = A.__cache;
    i.bindBuffer(i.UNIFORM_BUFFER, _);
    for (let L = 0, w = T.length; L < w; L++) {
      const V = T[L];
      if (m(V, L, C) === !0) {
        const M = V.__offset, y = Array.isArray(V.value) ? V.value : [V.value];
        let Y = 0;
        for (let ce = 0; ce < y.length; ce++) {
          const B = y[ce], H = x(B);
          typeof B == "number" ? (V.__data[0] = B, i.bufferSubData(i.UNIFORM_BUFFER, M + Y, V.__data)) : B.isMatrix3 ? (V.__data[0] = B.elements[0], V.__data[1] = B.elements[1], V.__data[2] = B.elements[2], V.__data[3] = B.elements[0], V.__data[4] = B.elements[3], V.__data[5] = B.elements[4], V.__data[6] = B.elements[5], V.__data[7] = B.elements[0], V.__data[8] = B.elements[6], V.__data[9] = B.elements[7], V.__data[10] = B.elements[8], V.__data[11] = B.elements[0]) : (B.toArray(V.__data, Y), Y += H.storage / Float32Array.BYTES_PER_ELEMENT);
        }
        i.bufferSubData(i.UNIFORM_BUFFER, M, V.__data);
      }
    }
    i.bindBuffer(i.UNIFORM_BUFFER, null);
  }
  function m(A, _, T) {
    const C = A.value;
    if (T[_] === void 0) {
      if (typeof C == "number")
        T[_] = C;
      else {
        const L = Array.isArray(C) ? C : [C], w = [];
        for (let V = 0; V < L.length; V++)
          w.push(L[V].clone());
        T[_] = w;
      }
      return !0;
    } else if (typeof C == "number") {
      if (T[_] !== C)
        return T[_] = C, !0;
    } else {
      const L = Array.isArray(T[_]) ? T[_] : [T[_]], w = Array.isArray(C) ? C : [C];
      for (let V = 0; V < L.length; V++) {
        const M = L[V];
        if (M.equals(w[V]) === !1)
          return M.copy(w[V]), !0;
      }
    }
    return !1;
  }
  function g(A) {
    const _ = A.uniforms;
    let T = 0;
    const C = 16;
    let L = 0;
    for (let w = 0, V = _.length; w < V; w++) {
      const M = _[w], y = {
        boundary: 0,
        // bytes
        storage: 0
        // bytes
      }, Y = Array.isArray(M.value) ? M.value : [M.value];
      for (let ce = 0, B = Y.length; ce < B; ce++) {
        const H = Y[ce], G = x(H);
        y.boundary += G.boundary, y.storage += G.storage;
      }
      if (M.__data = new Float32Array(y.storage / Float32Array.BYTES_PER_ELEMENT), M.__offset = T, w > 0) {
        L = T % C;
        const ce = C - L;
        L !== 0 && ce - y.boundary < 0 && (T += C - L, M.__offset = T);
      }
      T += y.storage;
    }
    return L = T % C, L > 0 && (T += C - L), A.__size = T, A.__cache = {}, this;
  }
  function x(A) {
    const _ = {
      boundary: 0,
      // bytes
      storage: 0
      // bytes
    };
    return typeof A == "number" ? (_.boundary = 4, _.storage = 4) : A.isVector2 ? (_.boundary = 8, _.storage = 8) : A.isVector3 || A.isColor ? (_.boundary = 16, _.storage = 12) : A.isVector4 ? (_.boundary = 16, _.storage = 16) : A.isMatrix3 ? (_.boundary = 48, _.storage = 48) : A.isMatrix4 ? (_.boundary = 64, _.storage = 64) : A.isTexture ? console.warn("THREE.WebGLRenderer: Texture samplers can not be part of an uniforms group.") : console.warn("THREE.WebGLRenderer: Unsupported uniform value type.", A), _;
  }
  function p(A) {
    const _ = A.target;
    _.removeEventListener("dispose", p);
    const T = o.indexOf(_.__bindingPointIndex);
    o.splice(T, 1), i.deleteBuffer(r[_.id]), delete r[_.id], delete s[_.id];
  }
  function d() {
    for (const A in r)
      i.deleteBuffer(r[A]);
    o = [], r = {}, s = {};
  }
  return {
    bind: l,
    update: c,
    dispose: d
  };
}
function dp() {
  const i = ir("canvas");
  return i.style.display = "block", i;
}
class Po {
  constructor(e = {}) {
    const {
      canvas: t = dp(),
      context: n = null,
      depth: r = !0,
      stencil: s = !0,
      alpha: o = !1,
      antialias: a = !1,
      premultipliedAlpha: l = !0,
      preserveDrawingBuffer: c = !1,
      powerPreference: h = "default",
      failIfMajorPerformanceCaveat: f = !1
    } = e;
    this.isWebGLRenderer = !0;
    let u;
    n !== null ? u = n.getContextAttributes().alpha : u = o;
    const m = new Uint32Array(4), g = new Int32Array(4);
    let x = null, p = null;
    const d = [], A = [];
    this.domElement = t, this.debug = {
      /**
       * Enables error checking and reporting when shader programs are being compiled
       * @type {boolean}
       */
      checkShaderErrors: !0,
      /**
       * Callback for custom error reporting.
       * @type {?Function}
       */
      onShaderError: null
    }, this.autoClear = !0, this.autoClearColor = !0, this.autoClearDepth = !0, this.autoClearStencil = !0, this.sortObjects = !0, this.clippingPlanes = [], this.localClippingEnabled = !1, this.outputColorSpace = Oe, this._useLegacyLights = !1, this.toneMapping = an, this.toneMappingExposure = 1;
    const _ = this;
    let T = !1, C = 0, L = 0, w = null, V = -1, M = null;
    const y = new ot(), Y = new ot();
    let ce = null;
    const B = new We(0);
    let H = 0, G = t.width, Q = t.height, X = 1, j = null, J = null;
    const ee = new ot(0, 0, G, Q), I = new ot(0, 0, G, Q);
    let q = !1;
    const pe = new us();
    let me = !1, ve = !1, Ae = null;
    const be = new nt(), Te = new oe(), ke = new U(), Ye = { background: null, fog: null, environment: null, overrideMaterial: null, isScene: !0 };
    function we() {
      return w === null ? X : 1;
    }
    let b = n;
    function le(S, N) {
      for (let W = 0; W < S.length; W++) {
        const D = S[W], k = t.getContext(D, N);
        if (k !== null)
          return k;
      }
      return null;
    }
    try {
      const S = {
        alpha: !0,
        depth: r,
        stencil: s,
        antialias: a,
        premultipliedAlpha: l,
        preserveDrawingBuffer: c,
        powerPreference: h,
        failIfMajorPerformanceCaveat: f
      };
      if ("setAttribute" in t && t.setAttribute("data-engine", `three.js r${ss}`), t.addEventListener("webglcontextlost", fe, !1), t.addEventListener("webglcontextrestored", F, !1), t.addEventListener("webglcontextcreationerror", ne, !1), b === null) {
        const N = ["webgl2", "webgl", "experimental-webgl"];
        if (_.isWebGL1Renderer === !0 && N.shift(), b = le(N, S), b === null)
          throw le(N) ? new Error("Error creating WebGL context with your selected attributes.") : new Error("Error creating WebGL context.");
      }
      typeof WebGLRenderingContext < "u" && b instanceof WebGLRenderingContext && console.warn("THREE.WebGLRenderer: WebGL 1 support was deprecated in r153 and will be removed in r163."), b.getShaderPrecisionFormat === void 0 && (b.getShaderPrecisionFormat = function() {
        return { rangeMin: 1, rangeMax: 1, precision: 1 };
      });
    } catch (S) {
      throw console.error("THREE.WebGLRenderer: " + S.message), S;
    }
    let Z, re, $, ye, xe, Me, Ce, ze, Ze, E, v, O, ie, te, se, Ee, ae, z, R, K, _e, he, ge, De;
    function Ve() {
      Z = new yf(b), re = new _f(b, Z, e), Z.init(re), he = new ap(b, Z, re), $ = new rp(b, Z, re), ye = new bf(b), xe = new Wd(), Me = new sp(b, Z, $, xe, re, he, ye), Ce = new xf(_), ze = new Ef(_), Ze = new Nc(b, re), ge = new mf(b, Z, Ze, re), E = new Tf(b, Ze, ye, ge), v = new Pf(b, E, Ze, ye), R = new Cf(b, re, Me), Ee = new vf(xe), O = new kd(_, Ce, ze, Z, re, ge, Ee), ie = new up(_, xe), te = new Yd(), se = new $d(Z, re), z = new pf(_, Ce, ze, $, v, u, l), ae = new ip(_, v, re), De = new fp(b, ye, re, $), K = new gf(b, Z, ye, re), _e = new Af(b, Z, ye, re), ye.programs = O.programs, _.capabilities = re, _.extensions = Z, _.properties = xe, _.renderLists = te, _.shadowMap = ae, _.state = $, _.info = ye;
    }
    Ve();
    const P = new hp(_, b);
    this.xr = P, this.getContext = function() {
      return b;
    }, this.getContextAttributes = function() {
      return b.getContextAttributes();
    }, this.forceContextLoss = function() {
      const S = Z.get("WEBGL_lose_context");
      S && S.loseContext();
    }, this.forceContextRestore = function() {
      const S = Z.get("WEBGL_lose_context");
      S && S.restoreContext();
    }, this.getPixelRatio = function() {
      return X;
    }, this.setPixelRatio = function(S) {
      S !== void 0 && (X = S, this.setSize(G, Q, !1));
    }, this.getSize = function(S) {
      return S.set(G, Q);
    }, this.setSize = function(S, N, W = !0) {
      if (P.isPresenting) {
        console.warn("THREE.WebGLRenderer: Can't change size while VR device is presenting.");
        return;
      }
      G = S, Q = N, t.width = Math.floor(S * X), t.height = Math.floor(N * X), W === !0 && (t.style.width = S + "px", t.style.height = N + "px"), this.setViewport(0, 0, S, N);
    }, this.getDrawingBufferSize = function(S) {
      return S.set(G * X, Q * X).floor();
    }, this.setDrawingBufferSize = function(S, N, W) {
      G = S, Q = N, X = W, t.width = Math.floor(S * W), t.height = Math.floor(N * W), this.setViewport(0, 0, S, N);
    }, this.getCurrentViewport = function(S) {
      return S.copy(y);
    }, this.getViewport = function(S) {
      return S.copy(ee);
    }, this.setViewport = function(S, N, W, D) {
      S.isVector4 ? ee.set(S.x, S.y, S.z, S.w) : ee.set(S, N, W, D), $.viewport(y.copy(ee).multiplyScalar(X).floor());
    }, this.getScissor = function(S) {
      return S.copy(I);
    }, this.setScissor = function(S, N, W, D) {
      S.isVector4 ? I.set(S.x, S.y, S.z, S.w) : I.set(S, N, W, D), $.scissor(Y.copy(I).multiplyScalar(X).floor());
    }, this.getScissorTest = function() {
      return q;
    }, this.setScissorTest = function(S) {
      $.setScissorTest(q = S);
    }, this.setOpaqueSort = function(S) {
      j = S;
    }, this.setTransparentSort = function(S) {
      J = S;
    }, this.getClearColor = function(S) {
      return S.copy(z.getClearColor());
    }, this.setClearColor = function() {
      z.setClearColor.apply(z, arguments);
    }, this.getClearAlpha = function() {
      return z.getClearAlpha();
    }, this.setClearAlpha = function() {
      z.setClearAlpha.apply(z, arguments);
    }, this.clear = function(S = !0, N = !0, W = !0) {
      let D = 0;
      if (S) {
        let k = !1;
        if (w !== null) {
          const Se = w.texture.format;
          k = Se === lo || Se === oo || Se === ao;
        }
        if (k) {
          const Se = w.texture.type, Re = Se === on || Se === nn || Se === as || Se === xn || Se === ro || Se === so, Le = z.getClearColor(), Ue = z.getClearAlpha(), He = Le.r, Pe = Le.g, Ie = Le.b;
          Re ? (m[0] = He, m[1] = Pe, m[2] = Ie, m[3] = Ue, b.clearBufferuiv(b.COLOR, 0, m)) : (g[0] = He, g[1] = Pe, g[2] = Ie, g[3] = Ue, b.clearBufferiv(b.COLOR, 0, g));
        } else
          D |= b.COLOR_BUFFER_BIT;
      }
      N && (D |= b.DEPTH_BUFFER_BIT), W && (D |= b.STENCIL_BUFFER_BIT), b.clear(D);
    }, this.clearColor = function() {
      this.clear(!0, !1, !1);
    }, this.clearDepth = function() {
      this.clear(!1, !0, !1);
    }, this.clearStencil = function() {
      this.clear(!1, !1, !0);
    }, this.dispose = function() {
      t.removeEventListener("webglcontextlost", fe, !1), t.removeEventListener("webglcontextrestored", F, !1), t.removeEventListener("webglcontextcreationerror", ne, !1), te.dispose(), se.dispose(), xe.dispose(), Ce.dispose(), ze.dispose(), v.dispose(), ge.dispose(), De.dispose(), O.dispose(), P.dispose(), P.removeEventListener("sessionstart", je), P.removeEventListener("sessionend", Dt), Ae && (Ae.dispose(), Ae = null), ut.stop();
    };
    function fe(S) {
      S.preventDefault(), console.log("THREE.WebGLRenderer: Context Lost."), T = !0;
    }
    function F() {
      console.log("THREE.WebGLRenderer: Context Restored."), T = !1;
      const S = ye.autoReset, N = ae.enabled, W = ae.autoUpdate, D = ae.needsUpdate, k = ae.type;
      Ve(), ye.autoReset = S, ae.enabled = N, ae.autoUpdate = W, ae.needsUpdate = D, ae.type = k;
    }
    function ne(S) {
      console.error("THREE.WebGLRenderer: A WebGL context could not be created. Reason: ", S.statusMessage);
    }
    function de(S) {
      const N = S.target;
      N.removeEventListener("dispose", de), Fe(N);
    }
    function Fe(S) {
      Xe(S), xe.remove(S);
    }
    function Xe(S) {
      const N = xe.get(S).programs;
      N !== void 0 && (N.forEach(function(W) {
        O.releaseProgram(W);
      }), S.isShaderMaterial && O.releaseShaderCache(S));
    }
    this.renderBufferDirect = function(S, N, W, D, k, Se) {
      N === null && (N = Ye);
      const Re = k.isMesh && k.matrixWorld.determinant() < 0, Le = Xo(S, N, W, D, k);
      $.setMaterial(D, Re);
      let Ue = W.index, He = 1;
      if (D.wireframe === !0) {
        if (Ue = E.getWireframeAttribute(W), Ue === void 0)
          return;
        He = 2;
      }
      const Pe = W.drawRange, Ie = W.attributes.position;
      let Je = Pe.start * He, $e = (Pe.start + Pe.count) * He;
      Se !== null && (Je = Math.max(Je, Se.start * He), $e = Math.min($e, (Se.start + Se.count) * He)), Ue !== null ? (Je = Math.max(Je, 0), $e = Math.min($e, Ue.count)) : Ie != null && (Je = Math.max(Je, 0), $e = Math.min($e, Ie.count));
      const Et = $e - Je;
      if (Et < 0 || Et === 1 / 0)
        return;
      ge.setup(k, D, Le, W, Ue);
      let zt, Qe = K;
      if (Ue !== null && (zt = Ze.get(Ue), Qe = _e, Qe.setIndex(zt)), k.isMesh)
        D.wireframe === !0 ? ($.setLineWidth(D.wireframeLinewidth * we()), Qe.setMode(b.LINES)) : Qe.setMode(b.TRIANGLES);
      else if (k.isLine) {
        let Ge = D.linewidth;
        Ge === void 0 && (Ge = 1), $.setLineWidth(Ge * we()), k.isLineSegments ? Qe.setMode(b.LINES) : k.isLineLoop ? Qe.setMode(b.LINE_LOOP) : Qe.setMode(b.LINE_STRIP);
      } else
        k.isPoints ? Qe.setMode(b.POINTS) : k.isSprite && Qe.setMode(b.TRIANGLES);
      if (k.isInstancedMesh)
        Qe.renderInstances(Je, Et, k.count);
      else if (W.isInstancedBufferGeometry) {
        const Ge = W._maxInstanceCount !== void 0 ? W._maxInstanceCount : 1 / 0, cr = Math.min(W.instanceCount, Ge);
        Qe.renderInstances(Je, Et, cr);
      } else
        Qe.render(Je, Et);
    }, this.compile = function(S, N) {
      function W(D, k, Se) {
        D.transparent === !0 && D.side === Mt && D.forceSinglePass === !1 ? (D.side = gt, D.needsUpdate = !0, Ri(D, k, Se), D.side = ln, D.needsUpdate = !0, Ri(D, k, Se), D.side = Mt) : Ri(D, k, Se);
      }
      p = se.get(S), p.init(), A.push(p), S.traverseVisible(function(D) {
        D.isLight && D.layers.test(N.layers) && (p.pushLight(D), D.castShadow && p.pushShadow(D));
      }), p.setupLights(_._useLegacyLights), S.traverse(function(D) {
        const k = D.material;
        if (k)
          if (Array.isArray(k))
            for (let Se = 0; Se < k.length; Se++) {
              const Re = k[Se];
              W(Re, S, D);
            }
          else
            W(k, S, D);
      }), A.pop(), p = null;
    };
    let qe = null;
    function Zt(S) {
      qe && qe(S);
    }
    function je() {
      ut.stop();
    }
    function Dt() {
      ut.start();
    }
    const ut = new To();
    ut.setAnimationLoop(Zt), typeof self < "u" && ut.setContext(self), this.setAnimationLoop = function(S) {
      qe = S, P.setAnimationLoop(S), S === null ? ut.stop() : ut.start();
    }, P.addEventListener("sessionstart", je), P.addEventListener("sessionend", Dt), this.render = function(S, N) {
      if (N !== void 0 && N.isCamera !== !0) {
        console.error("THREE.WebGLRenderer.render: camera is not an instance of THREE.Camera.");
        return;
      }
      if (T === !0)
        return;
      S.matrixWorldAutoUpdate === !0 && S.updateMatrixWorld(), N.parent === null && N.matrixWorldAutoUpdate === !0 && N.updateMatrixWorld(), P.enabled === !0 && P.isPresenting === !0 && (P.cameraAutoUpdate === !0 && P.updateCamera(N), N = P.getCamera()), S.isScene === !0 && S.onBeforeRender(_, S, N, w), p = se.get(S, A.length), p.init(), A.push(p), be.multiplyMatrices(N.projectionMatrix, N.matrixWorldInverse), pe.setFromProjectionMatrix(be), ve = this.localClippingEnabled, me = Ee.init(this.clippingPlanes, ve), x = te.get(S, d.length), x.init(), d.push(x), _s(S, N, 0, _.sortObjects), x.finish(), _.sortObjects === !0 && x.sort(j, J), this.info.render.frame++, me === !0 && Ee.beginShadows();
      const W = p.state.shadowsArray;
      if (ae.render(W, S, N), me === !0 && Ee.endShadows(), this.info.autoReset === !0 && this.info.reset(), z.render(x, S), p.setupLights(_._useLegacyLights), N.isArrayCamera) {
        const D = N.cameras;
        for (let k = 0, Se = D.length; k < Se; k++) {
          const Re = D[k];
          vs(x, S, Re, Re.viewport);
        }
      } else
        vs(x, S, N);
      w !== null && (Me.updateMultisampleRenderTarget(w), Me.updateRenderTargetMipmap(w)), S.isScene === !0 && S.onAfterRender(_, S, N), ge.resetDefaultState(), V = -1, M = null, A.pop(), A.length > 0 ? p = A[A.length - 1] : p = null, d.pop(), d.length > 0 ? x = d[d.length - 1] : x = null;
    };
    function _s(S, N, W, D) {
      if (S.visible === !1)
        return;
      if (S.layers.test(N.layers)) {
        if (S.isGroup)
          W = S.renderOrder;
        else if (S.isLOD)
          S.autoUpdate === !0 && S.update(N);
        else if (S.isLight)
          p.pushLight(S), S.castShadow && p.pushShadow(S);
        else if (S.isSprite) {
          if (!S.frustumCulled || pe.intersectsSprite(S)) {
            D && ke.setFromMatrixPosition(S.matrixWorld).applyMatrix4(be);
            const Re = v.update(S), Le = S.material;
            Le.visible && x.push(S, Re, Le, W, ke.z, null);
          }
        } else if ((S.isMesh || S.isLine || S.isPoints) && (!S.frustumCulled || pe.intersectsObject(S))) {
          const Re = v.update(S), Le = S.material;
          if (D && (S.boundingSphere !== void 0 ? (S.boundingSphere === null && S.computeBoundingSphere(), ke.copy(S.boundingSphere.center)) : (Re.boundingSphere === null && Re.computeBoundingSphere(), ke.copy(Re.boundingSphere.center)), ke.applyMatrix4(S.matrixWorld).applyMatrix4(be)), Array.isArray(Le)) {
            const Ue = Re.groups;
            for (let He = 0, Pe = Ue.length; He < Pe; He++) {
              const Ie = Ue[He], Je = Le[Ie.materialIndex];
              Je && Je.visible && x.push(S, Re, Je, W, ke.z, Ie);
            }
          } else
            Le.visible && x.push(S, Re, Le, W, ke.z, null);
        }
      }
      const Se = S.children;
      for (let Re = 0, Le = Se.length; Re < Le; Re++)
        _s(Se[Re], N, W, D);
    }
    function vs(S, N, W, D) {
      const k = S.opaque, Se = S.transmissive, Re = S.transparent;
      p.setupLightsView(W), me === !0 && Ee.setGlobalState(_.clippingPlanes, W), Se.length > 0 && Wo(k, Se, N, W), D && $.viewport(y.copy(D)), k.length > 0 && wi(k, N, W), Se.length > 0 && wi(Se, N, W), Re.length > 0 && wi(Re, N, W), $.buffers.depth.setTest(!0), $.buffers.depth.setMask(!0), $.buffers.color.setMask(!0), $.setPolygonOffset(!1);
    }
    function Wo(S, N, W, D) {
      const k = re.isWebGL2;
      Ae === null && (Ae = new Tn(1, 1, {
        generateMipmaps: !0,
        type: Z.has("EXT_color_buffer_half_float") ? xi : on,
        minFilter: vi,
        samples: k ? 4 : 0
      })), _.getDrawingBufferSize(Te), k ? Ae.setSize(Te.x, Te.y) : Ae.setSize(nr(Te.x), nr(Te.y));
      const Se = _.getRenderTarget();
      _.setRenderTarget(Ae), _.getClearColor(B), H = _.getClearAlpha(), H < 1 && _.setClearColor(16777215, 0.5), _.clear();
      const Re = _.toneMapping;
      _.toneMapping = an, wi(S, W, D), Me.updateMultisampleRenderTarget(Ae), Me.updateRenderTargetMipmap(Ae);
      let Le = !1;
      for (let Ue = 0, He = N.length; Ue < He; Ue++) {
        const Pe = N[Ue], Ie = Pe.object, Je = Pe.geometry, $e = Pe.material, Et = Pe.group;
        if ($e.side === Mt && Ie.layers.test(D.layers)) {
          const zt = $e.side;
          $e.side = gt, $e.needsUpdate = !0, xs(Ie, W, D, Je, $e, Et), $e.side = zt, $e.needsUpdate = !0, Le = !0;
        }
      }
      Le === !0 && (Me.updateMultisampleRenderTarget(Ae), Me.updateRenderTargetMipmap(Ae)), _.setRenderTarget(Se), _.setClearColor(B, H), _.toneMapping = Re;
    }
    function wi(S, N, W) {
      const D = N.isScene === !0 ? N.overrideMaterial : null;
      for (let k = 0, Se = S.length; k < Se; k++) {
        const Re = S[k], Le = Re.object, Ue = Re.geometry, He = D === null ? Re.material : D, Pe = Re.group;
        Le.layers.test(W.layers) && xs(Le, N, W, Ue, He, Pe);
      }
    }
    function xs(S, N, W, D, k, Se) {
      S.onBeforeRender(_, N, W, D, k, Se), S.modelViewMatrix.multiplyMatrices(W.matrixWorldInverse, S.matrixWorld), S.normalMatrix.getNormalMatrix(S.modelViewMatrix), k.onBeforeRender(_, N, W, D, S, Se), k.transparent === !0 && k.side === Mt && k.forceSinglePass === !1 ? (k.side = gt, k.needsUpdate = !0, _.renderBufferDirect(W, N, D, k, S, Se), k.side = ln, k.needsUpdate = !0, _.renderBufferDirect(W, N, D, k, S, Se), k.side = Mt) : _.renderBufferDirect(W, N, D, k, S, Se), S.onAfterRender(_, N, W, D, k, Se);
    }
    function Ri(S, N, W) {
      N.isScene !== !0 && (N = Ye);
      const D = xe.get(S), k = p.state.lights, Se = p.state.shadowsArray, Re = k.state.version, Le = O.getParameters(S, k.state, Se, N, W), Ue = O.getProgramCacheKey(Le);
      let He = D.programs;
      D.environment = S.isMeshStandardMaterial ? N.environment : null, D.fog = N.fog, D.envMap = (S.isMeshStandardMaterial ? ze : Ce).get(S.envMap || D.environment), He === void 0 && (S.addEventListener("dispose", de), He = /* @__PURE__ */ new Map(), D.programs = He);
      let Pe = He.get(Ue);
      if (Pe !== void 0) {
        if (D.currentProgram === Pe && D.lightsStateVersion === Re)
          return Ms(S, Le), Pe;
      } else
        Le.uniforms = O.getUniforms(S), S.onBuild(W, Le, _), S.onBeforeCompile(Le, _), Pe = O.acquireProgram(Le, Ue), He.set(Ue, Pe), D.uniforms = Le.uniforms;
      const Ie = D.uniforms;
      (!S.isShaderMaterial && !S.isRawShaderMaterial || S.clipping === !0) && (Ie.clippingPlanes = Ee.uniform), Ms(S, Le), D.needsLights = qo(S), D.lightsStateVersion = Re, D.needsLights && (Ie.ambientLightColor.value = k.state.ambient, Ie.lightProbe.value = k.state.probe, Ie.directionalLights.value = k.state.directional, Ie.directionalLightShadows.value = k.state.directionalShadow, Ie.spotLights.value = k.state.spot, Ie.spotLightShadows.value = k.state.spotShadow, Ie.rectAreaLights.value = k.state.rectArea, Ie.ltc_1.value = k.state.rectAreaLTC1, Ie.ltc_2.value = k.state.rectAreaLTC2, Ie.pointLights.value = k.state.point, Ie.pointLightShadows.value = k.state.pointShadow, Ie.hemisphereLights.value = k.state.hemi, Ie.directionalShadowMap.value = k.state.directionalShadowMap, Ie.directionalShadowMatrix.value = k.state.directionalShadowMatrix, Ie.spotShadowMap.value = k.state.spotShadowMap, Ie.spotLightMatrix.value = k.state.spotLightMatrix, Ie.spotLightMap.value = k.state.spotLightMap, Ie.pointShadowMap.value = k.state.pointShadowMap, Ie.pointShadowMatrix.value = k.state.pointShadowMatrix);
      const Je = Pe.getUniforms(), $e = er.seqWithValue(Je.seq, Ie);
      return D.currentProgram = Pe, D.uniformsList = $e, Pe;
    }
    function Ms(S, N) {
      const W = xe.get(S);
      W.outputColorSpace = N.outputColorSpace, W.instancing = N.instancing, W.instancingColor = N.instancingColor, W.skinning = N.skinning, W.morphTargets = N.morphTargets, W.morphNormals = N.morphNormals, W.morphColors = N.morphColors, W.morphTargetsCount = N.morphTargetsCount, W.numClippingPlanes = N.numClippingPlanes, W.numIntersection = N.numClipIntersection, W.vertexAlphas = N.vertexAlphas, W.vertexTangents = N.vertexTangents, W.toneMapping = N.toneMapping;
    }
    function Xo(S, N, W, D, k) {
      N.isScene !== !0 && (N = Ye), Me.resetTextureUnits();
      const Se = N.fog, Re = D.isMeshStandardMaterial ? N.environment : null, Le = w === null ? _.outputColorSpace : w.isXRRenderTarget === !0 ? w.texture.colorSpace : Ft, Ue = (D.isMeshStandardMaterial ? ze : Ce).get(D.envMap || Re), He = D.vertexColors === !0 && !!W.attributes.color && W.attributes.color.itemSize === 4, Pe = !!W.attributes.tangent && (!!D.normalMap || D.anisotropy > 0), Ie = !!W.morphAttributes.position, Je = !!W.morphAttributes.normal, $e = !!W.morphAttributes.color;
      let Et = an;
      D.toneMapped && (w === null || w.isXRRenderTarget === !0) && (Et = _.toneMapping);
      const zt = W.morphAttributes.position || W.morphAttributes.normal || W.morphAttributes.color, Qe = zt !== void 0 ? zt.length : 0, Ge = xe.get(D), cr = p.state.lights;
      if (me === !0 && (ve === !0 || S !== M)) {
        const _t = S === M && D.id === V;
        Ee.setState(D, S, _t);
      }
      let et = !1;
      D.version === Ge.__version ? (Ge.needsLights && Ge.lightsStateVersion !== cr.state.version || Ge.outputColorSpace !== Le || k.isInstancedMesh && Ge.instancing === !1 || !k.isInstancedMesh && Ge.instancing === !0 || k.isSkinnedMesh && Ge.skinning === !1 || !k.isSkinnedMesh && Ge.skinning === !0 || k.isInstancedMesh && Ge.instancingColor === !0 && k.instanceColor === null || k.isInstancedMesh && Ge.instancingColor === !1 && k.instanceColor !== null || Ge.envMap !== Ue || D.fog === !0 && Ge.fog !== Se || Ge.numClippingPlanes !== void 0 && (Ge.numClippingPlanes !== Ee.numPlanes || Ge.numIntersection !== Ee.numIntersection) || Ge.vertexAlphas !== He || Ge.vertexTangents !== Pe || Ge.morphTargets !== Ie || Ge.morphNormals !== Je || Ge.morphColors !== $e || Ge.toneMapping !== Et || re.isWebGL2 === !0 && Ge.morphTargetsCount !== Qe) && (et = !0) : (et = !0, Ge.__version = D.version);
      let hn = Ge.currentProgram;
      et === !0 && (hn = Ri(D, N, k));
      let Ss = !1, si = !1, hr = !1;
      const ft = hn.getUniforms(), un = Ge.uniforms;
      if ($.useProgram(hn.program) && (Ss = !0, si = !0, hr = !0), D.id !== V && (V = D.id, si = !0), Ss || M !== S) {
        if (ft.setValue(b, "projectionMatrix", S.projectionMatrix), re.logarithmicDepthBuffer && ft.setValue(
          b,
          "logDepthBufFC",
          2 / (Math.log(S.far + 1) / Math.LN2)
        ), M !== S && (M = S, si = !0, hr = !0), D.isShaderMaterial || D.isMeshPhongMaterial || D.isMeshToonMaterial || D.isMeshStandardMaterial || D.envMap) {
          const _t = ft.map.cameraPosition;
          _t !== void 0 && _t.setValue(
            b,
            ke.setFromMatrixPosition(S.matrixWorld)
          );
        }
        (D.isMeshPhongMaterial || D.isMeshToonMaterial || D.isMeshLambertMaterial || D.isMeshBasicMaterial || D.isMeshStandardMaterial || D.isShaderMaterial) && ft.setValue(b, "isOrthographic", S.isOrthographicCamera === !0), (D.isMeshPhongMaterial || D.isMeshToonMaterial || D.isMeshLambertMaterial || D.isMeshBasicMaterial || D.isMeshStandardMaterial || D.isShaderMaterial || D.isShadowMaterial || k.isSkinnedMesh) && ft.setValue(b, "viewMatrix", S.matrixWorldInverse);
      }
      if (k.isSkinnedMesh) {
        ft.setOptional(b, k, "bindMatrix"), ft.setOptional(b, k, "bindMatrixInverse");
        const _t = k.skeleton;
        _t && (re.floatVertexTextures ? (_t.boneTexture === null && _t.computeBoneTexture(), ft.setValue(b, "boneTexture", _t.boneTexture, Me), ft.setValue(b, "boneTextureSize", _t.boneTextureSize)) : console.warn("THREE.WebGLRenderer: SkinnedMesh can only be used with WebGL 2. With WebGL 1 OES_texture_float and vertex textures support is required."));
      }
      const ur = W.morphAttributes;
      if ((ur.position !== void 0 || ur.normal !== void 0 || ur.color !== void 0 && re.isWebGL2 === !0) && R.update(k, W, hn), (si || Ge.receiveShadow !== k.receiveShadow) && (Ge.receiveShadow = k.receiveShadow, ft.setValue(b, "receiveShadow", k.receiveShadow)), D.isMeshGouraudMaterial && D.envMap !== null && (un.envMap.value = Ue, un.flipEnvMap.value = Ue.isCubeTexture && Ue.isRenderTargetTexture === !1 ? -1 : 1), si && (ft.setValue(b, "toneMappingExposure", _.toneMappingExposure), Ge.needsLights && Yo(un, hr), Se && D.fog === !0 && ie.refreshFogUniforms(un, Se), ie.refreshMaterialUniforms(un, D, X, Q, Ae), er.upload(b, Ge.uniformsList, un, Me)), D.isShaderMaterial && D.uniformsNeedUpdate === !0 && (er.upload(b, Ge.uniformsList, un, Me), D.uniformsNeedUpdate = !1), D.isSpriteMaterial && ft.setValue(b, "center", k.center), ft.setValue(b, "modelViewMatrix", k.modelViewMatrix), ft.setValue(b, "normalMatrix", k.normalMatrix), ft.setValue(b, "modelMatrix", k.matrixWorld), D.isShaderMaterial || D.isRawShaderMaterial) {
        const _t = D.uniformsGroups;
        for (let fr = 0, jo = _t.length; fr < jo; fr++)
          if (re.isWebGL2) {
            const Es = _t[fr];
            De.update(Es, hn), De.bind(Es, hn);
          } else
            console.warn("THREE.WebGLRenderer: Uniform Buffer Objects can only be used with WebGL 2.");
      }
      return hn;
    }
    function Yo(S, N) {
      S.ambientLightColor.needsUpdate = N, S.lightProbe.needsUpdate = N, S.directionalLights.needsUpdate = N, S.directionalLightShadows.needsUpdate = N, S.pointLights.needsUpdate = N, S.pointLightShadows.needsUpdate = N, S.spotLights.needsUpdate = N, S.spotLightShadows.needsUpdate = N, S.rectAreaLights.needsUpdate = N, S.hemisphereLights.needsUpdate = N;
    }
    function qo(S) {
      return S.isMeshLambertMaterial || S.isMeshToonMaterial || S.isMeshPhongMaterial || S.isMeshStandardMaterial || S.isShadowMaterial || S.isShaderMaterial && S.lights === !0;
    }
    this.getActiveCubeFace = function() {
      return C;
    }, this.getActiveMipmapLevel = function() {
      return L;
    }, this.getRenderTarget = function() {
      return w;
    }, this.setRenderTargetTextures = function(S, N, W) {
      xe.get(S.texture).__webglTexture = N, xe.get(S.depthTexture).__webglTexture = W;
      const D = xe.get(S);
      D.__hasExternalTextures = !0, D.__hasExternalTextures && (D.__autoAllocateDepthBuffer = W === void 0, D.__autoAllocateDepthBuffer || Z.has("WEBGL_multisampled_render_to_texture") === !0 && (console.warn("THREE.WebGLRenderer: Render-to-texture extension was disabled because an external texture was provided"), D.__useRenderToTexture = !1));
    }, this.setRenderTargetFramebuffer = function(S, N) {
      const W = xe.get(S);
      W.__webglFramebuffer = N, W.__useDefaultFramebuffer = N === void 0;
    }, this.setRenderTarget = function(S, N = 0, W = 0) {
      w = S, C = N, L = W;
      let D = !0, k = null, Se = !1, Re = !1;
      if (S) {
        const Ue = xe.get(S);
        Ue.__useDefaultFramebuffer !== void 0 ? ($.bindFramebuffer(b.FRAMEBUFFER, null), D = !1) : Ue.__webglFramebuffer === void 0 ? Me.setupRenderTarget(S) : Ue.__hasExternalTextures && Me.rebindTextures(S, xe.get(S.texture).__webglTexture, xe.get(S.depthTexture).__webglTexture);
        const He = S.texture;
        (He.isData3DTexture || He.isDataArrayTexture || He.isCompressedArrayTexture) && (Re = !0);
        const Pe = xe.get(S).__webglFramebuffer;
        S.isWebGLCubeRenderTarget ? (Array.isArray(Pe[N]) ? k = Pe[N][W] : k = Pe[N], Se = !0) : re.isWebGL2 && S.samples > 0 && Me.useMultisampledRTT(S) === !1 ? k = xe.get(S).__webglMultisampledFramebuffer : Array.isArray(Pe) ? k = Pe[W] : k = Pe, y.copy(S.viewport), Y.copy(S.scissor), ce = S.scissorTest;
      } else
        y.copy(ee).multiplyScalar(X).floor(), Y.copy(I).multiplyScalar(X).floor(), ce = q;
      if ($.bindFramebuffer(b.FRAMEBUFFER, k) && re.drawBuffers && D && $.drawBuffers(S, k), $.viewport(y), $.scissor(Y), $.setScissorTest(ce), Se) {
        const Ue = xe.get(S.texture);
        b.framebufferTexture2D(b.FRAMEBUFFER, b.COLOR_ATTACHMENT0, b.TEXTURE_CUBE_MAP_POSITIVE_X + N, Ue.__webglTexture, W);
      } else if (Re) {
        const Ue = xe.get(S.texture), He = N || 0;
        b.framebufferTextureLayer(b.FRAMEBUFFER, b.COLOR_ATTACHMENT0, Ue.__webglTexture, W || 0, He);
      }
      V = -1;
    }, this.readRenderTargetPixels = function(S, N, W, D, k, Se, Re) {
      if (!(S && S.isWebGLRenderTarget)) {
        console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not THREE.WebGLRenderTarget.");
        return;
      }
      let Le = xe.get(S).__webglFramebuffer;
      if (S.isWebGLCubeRenderTarget && Re !== void 0 && (Le = Le[Re]), Le) {
        $.bindFramebuffer(b.FRAMEBUFFER, Le);
        try {
          const Ue = S.texture, He = Ue.format, Pe = Ue.type;
          if (He !== Ut && he.convert(He) !== b.getParameter(b.IMPLEMENTATION_COLOR_READ_FORMAT)) {
            console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not in RGBA or implementation defined format.");
            return;
          }
          const Ie = Pe === xi && (Z.has("EXT_color_buffer_half_float") || re.isWebGL2 && Z.has("EXT_color_buffer_float"));
          if (Pe !== on && he.convert(Pe) !== b.getParameter(b.IMPLEMENTATION_COLOR_READ_TYPE) && // Edge and Chrome Mac < 52 (#9513)
          !(Pe === rn && (re.isWebGL2 || Z.has("OES_texture_float") || Z.has("WEBGL_color_buffer_float"))) && // Chrome Mac >= 52 and Firefox
          !Ie) {
            console.error("THREE.WebGLRenderer.readRenderTargetPixels: renderTarget is not in UnsignedByteType or implementation defined type.");
            return;
          }
          N >= 0 && N <= S.width - D && W >= 0 && W <= S.height - k && b.readPixels(N, W, D, k, he.convert(He), he.convert(Pe), Se);
        } finally {
          const Ue = w !== null ? xe.get(w).__webglFramebuffer : null;
          $.bindFramebuffer(b.FRAMEBUFFER, Ue);
        }
      }
    }, this.copyFramebufferToTexture = function(S, N, W = 0) {
      const D = Math.pow(2, -W), k = Math.floor(N.image.width * D), Se = Math.floor(N.image.height * D);
      Me.setTexture2D(N, 0), b.copyTexSubImage2D(b.TEXTURE_2D, W, 0, 0, S.x, S.y, k, Se), $.unbindTexture();
    }, this.copyTextureToTexture = function(S, N, W, D = 0) {
      const k = N.image.width, Se = N.image.height, Re = he.convert(W.format), Le = he.convert(W.type);
      Me.setTexture2D(W, 0), b.pixelStorei(b.UNPACK_FLIP_Y_WEBGL, W.flipY), b.pixelStorei(b.UNPACK_PREMULTIPLY_ALPHA_WEBGL, W.premultiplyAlpha), b.pixelStorei(b.UNPACK_ALIGNMENT, W.unpackAlignment), N.isDataTexture ? b.texSubImage2D(b.TEXTURE_2D, D, S.x, S.y, k, Se, Re, Le, N.image.data) : N.isCompressedTexture ? b.compressedTexSubImage2D(b.TEXTURE_2D, D, S.x, S.y, N.mipmaps[0].width, N.mipmaps[0].height, Re, N.mipmaps[0].data) : b.texSubImage2D(b.TEXTURE_2D, D, S.x, S.y, Re, Le, N.image), D === 0 && W.generateMipmaps && b.generateMipmap(b.TEXTURE_2D), $.unbindTexture();
    }, this.copyTextureToTexture3D = function(S, N, W, D, k = 0) {
      if (_.isWebGL1Renderer) {
        console.warn("THREE.WebGLRenderer.copyTextureToTexture3D: can only be used with WebGL2.");
        return;
      }
      const Se = S.max.x - S.min.x + 1, Re = S.max.y - S.min.y + 1, Le = S.max.z - S.min.z + 1, Ue = he.convert(D.format), He = he.convert(D.type);
      let Pe;
      if (D.isData3DTexture)
        Me.setTexture3D(D, 0), Pe = b.TEXTURE_3D;
      else if (D.isDataArrayTexture)
        Me.setTexture2DArray(D, 0), Pe = b.TEXTURE_2D_ARRAY;
      else {
        console.warn("THREE.WebGLRenderer.copyTextureToTexture3D: only supports THREE.DataTexture3D and THREE.DataTexture2DArray.");
        return;
      }
      b.pixelStorei(b.UNPACK_FLIP_Y_WEBGL, D.flipY), b.pixelStorei(b.UNPACK_PREMULTIPLY_ALPHA_WEBGL, D.premultiplyAlpha), b.pixelStorei(b.UNPACK_ALIGNMENT, D.unpackAlignment);
      const Ie = b.getParameter(b.UNPACK_ROW_LENGTH), Je = b.getParameter(b.UNPACK_IMAGE_HEIGHT), $e = b.getParameter(b.UNPACK_SKIP_PIXELS), Et = b.getParameter(b.UNPACK_SKIP_ROWS), zt = b.getParameter(b.UNPACK_SKIP_IMAGES), Qe = W.isCompressedTexture ? W.mipmaps[0] : W.image;
      b.pixelStorei(b.UNPACK_ROW_LENGTH, Qe.width), b.pixelStorei(b.UNPACK_IMAGE_HEIGHT, Qe.height), b.pixelStorei(b.UNPACK_SKIP_PIXELS, S.min.x), b.pixelStorei(b.UNPACK_SKIP_ROWS, S.min.y), b.pixelStorei(b.UNPACK_SKIP_IMAGES, S.min.z), W.isDataTexture || W.isData3DTexture ? b.texSubImage3D(Pe, k, N.x, N.y, N.z, Se, Re, Le, Ue, He, Qe.data) : W.isCompressedArrayTexture ? (console.warn("THREE.WebGLRenderer.copyTextureToTexture3D: untested support for compressed srcTexture."), b.compressedTexSubImage3D(Pe, k, N.x, N.y, N.z, Se, Re, Le, Ue, Qe.data)) : b.texSubImage3D(Pe, k, N.x, N.y, N.z, Se, Re, Le, Ue, He, Qe), b.pixelStorei(b.UNPACK_ROW_LENGTH, Ie), b.pixelStorei(b.UNPACK_IMAGE_HEIGHT, Je), b.pixelStorei(b.UNPACK_SKIP_PIXELS, $e), b.pixelStorei(b.UNPACK_SKIP_ROWS, Et), b.pixelStorei(b.UNPACK_SKIP_IMAGES, zt), k === 0 && D.generateMipmaps && b.generateMipmap(Pe), $.unbindTexture();
    }, this.initTexture = function(S) {
      S.isCubeTexture ? Me.setTextureCube(S, 0) : S.isData3DTexture ? Me.setTexture3D(S, 0) : S.isDataArrayTexture || S.isCompressedArrayTexture ? Me.setTexture2DArray(S, 0) : Me.setTexture2D(S, 0), $.unbindTexture();
    }, this.resetState = function() {
      C = 0, L = 0, w = null, $.reset(), ge.reset();
    }, typeof __THREE_DEVTOOLS__ < "u" && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("observe", { detail: this }));
  }
  get coordinateSystem() {
    return qt;
  }
  get physicallyCorrectLights() {
    return console.warn("THREE.WebGLRenderer: The property .physicallyCorrectLights has been removed. Set renderer.useLegacyLights instead."), !this.useLegacyLights;
  }
  set physicallyCorrectLights(e) {
    console.warn("THREE.WebGLRenderer: The property .physicallyCorrectLights has been removed. Set renderer.useLegacyLights instead."), this.useLegacyLights = !e;
  }
  get outputEncoding() {
    return console.warn("THREE.WebGLRenderer: Property .outputEncoding has been removed. Use .outputColorSpace instead."), this.outputColorSpace === Oe ? Sn : co;
  }
  set outputEncoding(e) {
    console.warn("THREE.WebGLRenderer: Property .outputEncoding has been removed. Use .outputColorSpace instead."), this.outputColorSpace = e === Sn ? Oe : Ft;
  }
  get useLegacyLights() {
    return console.warn("THREE.WebGLRenderer: The property .useLegacyLights has been deprecated. Migrate your lighting according to the following guide: https://discourse.threejs.org/t/updates-to-lighting-in-three-js-r155/53733."), this._useLegacyLights;
  }
  set useLegacyLights(e) {
    console.warn("THREE.WebGLRenderer: The property .useLegacyLights has been deprecated. Migrate your lighting according to the following guide: https://discourse.threejs.org/t/updates-to-lighting-in-three-js-r155/53733."), this._useLegacyLights = e;
  }
}
class pp extends Po {
}
pp.prototype.isWebGL1Renderer = !0;
class mp extends ht {
  constructor() {
    super(), this.isScene = !0, this.type = "Scene", this.background = null, this.environment = null, this.fog = null, this.backgroundBlurriness = 0, this.backgroundIntensity = 1, this.overrideMaterial = null, typeof __THREE_DEVTOOLS__ < "u" && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("observe", { detail: this }));
  }
  copy(e, t) {
    return super.copy(e, t), e.background !== null && (this.background = e.background.clone()), e.environment !== null && (this.environment = e.environment.clone()), e.fog !== null && (this.fog = e.fog.clone()), this.backgroundBlurriness = e.backgroundBlurriness, this.backgroundIntensity = e.backgroundIntensity, e.overrideMaterial !== null && (this.overrideMaterial = e.overrideMaterial.clone()), this.matrixAutoUpdate = e.matrixAutoUpdate, this;
  }
  toJSON(e) {
    const t = super.toJSON(e);
    return this.fog !== null && (t.object.fog = this.fog.toJSON()), this.backgroundBlurriness > 0 && (t.object.backgroundBlurriness = this.backgroundBlurriness), this.backgroundIntensity !== 1 && (t.object.backgroundIntensity = this.backgroundIntensity), t;
  }
}
class Bt {
  constructor() {
    this.type = "Curve", this.arcLengthDivisions = 200;
  }
  // Virtual base class method to overwrite and implement in subclasses
  //	- t [0 .. 1]
  getPoint() {
    return console.warn("THREE.Curve: .getPoint() not implemented."), null;
  }
  // Get point at relative position in curve according to arc length
  // - u [0 .. 1]
  getPointAt(e, t) {
    const n = this.getUtoTmapping(e);
    return this.getPoint(n, t);
  }
  // Get sequence of points using getPoint( t )
  getPoints(e = 5) {
    const t = [];
    for (let n = 0; n <= e; n++)
      t.push(this.getPoint(n / e));
    return t;
  }
  // Get sequence of points using getPointAt( u )
  getSpacedPoints(e = 5) {
    const t = [];
    for (let n = 0; n <= e; n++)
      t.push(this.getPointAt(n / e));
    return t;
  }
  // Get total curve arc length
  getLength() {
    const e = this.getLengths();
    return e[e.length - 1];
  }
  // Get list of cumulative segment lengths
  getLengths(e = this.arcLengthDivisions) {
    if (this.cacheArcLengths && this.cacheArcLengths.length === e + 1 && !this.needsUpdate)
      return this.cacheArcLengths;
    this.needsUpdate = !1;
    const t = [];
    let n, r = this.getPoint(0), s = 0;
    t.push(0);
    for (let o = 1; o <= e; o++)
      n = this.getPoint(o / e), s += n.distanceTo(r), t.push(s), r = n;
    return this.cacheArcLengths = t, t;
  }
  updateArcLengths() {
    this.needsUpdate = !0, this.getLengths();
  }
  // Given u ( 0 .. 1 ), get a t to find p. This gives you points which are equidistant
  getUtoTmapping(e, t) {
    const n = this.getLengths();
    let r = 0;
    const s = n.length;
    let o;
    t ? o = t : o = e * n[s - 1];
    let a = 0, l = s - 1, c;
    for (; a <= l; )
      if (r = Math.floor(a + (l - a) / 2), c = n[r] - o, c < 0)
        a = r + 1;
      else if (c > 0)
        l = r - 1;
      else {
        l = r;
        break;
      }
    if (r = l, n[r] === o)
      return r / (s - 1);
    const h = n[r], u = n[r + 1] - h, m = (o - h) / u;
    return (r + m) / (s - 1);
  }
  // Returns a unit vector tangent at t
  // In case any sub curve does not implement its tangent derivation,
  // 2 points a small delta apart will be used to find its gradient
  // which seems to give a reasonable approximation
  getTangent(e, t) {
    let r = e - 1e-4, s = e + 1e-4;
    r < 0 && (r = 0), s > 1 && (s = 1);
    const o = this.getPoint(r), a = this.getPoint(s), l = t || (o.isVector2 ? new oe() : new U());
    return l.copy(a).sub(o).normalize(), l;
  }
  getTangentAt(e, t) {
    const n = this.getUtoTmapping(e);
    return this.getTangent(n, t);
  }
  computeFrenetFrames(e, t) {
    const n = new U(), r = [], s = [], o = [], a = new U(), l = new nt();
    for (let m = 0; m <= e; m++) {
      const g = m / e;
      r[m] = this.getTangentAt(g, new U());
    }
    s[0] = new U(), o[0] = new U();
    let c = Number.MAX_VALUE;
    const h = Math.abs(r[0].x), f = Math.abs(r[0].y), u = Math.abs(r[0].z);
    h <= c && (c = h, n.set(1, 0, 0)), f <= c && (c = f, n.set(0, 1, 0)), u <= c && n.set(0, 0, 1), a.crossVectors(r[0], n).normalize(), s[0].crossVectors(r[0], a), o[0].crossVectors(r[0], s[0]);
    for (let m = 1; m <= e; m++) {
      if (s[m] = s[m - 1].clone(), o[m] = o[m - 1].clone(), a.crossVectors(r[m - 1], r[m]), a.length() > Number.EPSILON) {
        a.normalize();
        const g = Math.acos(it(r[m - 1].dot(r[m]), -1, 1));
        s[m].applyMatrix4(l.makeRotationAxis(a, g));
      }
      o[m].crossVectors(r[m], s[m]);
    }
    if (t === !0) {
      let m = Math.acos(it(s[0].dot(s[e]), -1, 1));
      m /= e, r[0].dot(a.crossVectors(s[0], s[e])) > 0 && (m = -m);
      for (let g = 1; g <= e; g++)
        s[g].applyMatrix4(l.makeRotationAxis(r[g], m * g)), o[g].crossVectors(r[g], s[g]);
    }
    return {
      tangents: r,
      normals: s,
      binormals: o
    };
  }
  clone() {
    return new this.constructor().copy(this);
  }
  copy(e) {
    return this.arcLengthDivisions = e.arcLengthDivisions, this;
  }
  toJSON() {
    const e = {
      metadata: {
        version: 4.6,
        type: "Curve",
        generator: "Curve.toJSON"
      }
    };
    return e.arcLengthDivisions = this.arcLengthDivisions, e.type = this.type, e;
  }
  fromJSON(e) {
    return this.arcLengthDivisions = e.arcLengthDivisions, this;
  }
}
class ds extends Bt {
  constructor(e = 0, t = 0, n = 1, r = 1, s = 0, o = Math.PI * 2, a = !1, l = 0) {
    super(), this.isEllipseCurve = !0, this.type = "EllipseCurve", this.aX = e, this.aY = t, this.xRadius = n, this.yRadius = r, this.aStartAngle = s, this.aEndAngle = o, this.aClockwise = a, this.aRotation = l;
  }
  getPoint(e, t) {
    const n = t || new oe(), r = Math.PI * 2;
    let s = this.aEndAngle - this.aStartAngle;
    const o = Math.abs(s) < Number.EPSILON;
    for (; s < 0; )
      s += r;
    for (; s > r; )
      s -= r;
    s < Number.EPSILON && (o ? s = 0 : s = r), this.aClockwise === !0 && !o && (s === r ? s = -r : s = s - r);
    const a = this.aStartAngle + e * s;
    let l = this.aX + this.xRadius * Math.cos(a), c = this.aY + this.yRadius * Math.sin(a);
    if (this.aRotation !== 0) {
      const h = Math.cos(this.aRotation), f = Math.sin(this.aRotation), u = l - this.aX, m = c - this.aY;
      l = u * h - m * f + this.aX, c = u * f + m * h + this.aY;
    }
    return n.set(l, c);
  }
  copy(e) {
    return super.copy(e), this.aX = e.aX, this.aY = e.aY, this.xRadius = e.xRadius, this.yRadius = e.yRadius, this.aStartAngle = e.aStartAngle, this.aEndAngle = e.aEndAngle, this.aClockwise = e.aClockwise, this.aRotation = e.aRotation, this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.aX = this.aX, e.aY = this.aY, e.xRadius = this.xRadius, e.yRadius = this.yRadius, e.aStartAngle = this.aStartAngle, e.aEndAngle = this.aEndAngle, e.aClockwise = this.aClockwise, e.aRotation = this.aRotation, e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.aX = e.aX, this.aY = e.aY, this.xRadius = e.xRadius, this.yRadius = e.yRadius, this.aStartAngle = e.aStartAngle, this.aEndAngle = e.aEndAngle, this.aClockwise = e.aClockwise, this.aRotation = e.aRotation, this;
  }
}
class gp extends ds {
  constructor(e, t, n, r, s, o) {
    super(e, t, n, n, r, s, o), this.isArcCurve = !0, this.type = "ArcCurve";
  }
}
function ps() {
  let i = 0, e = 0, t = 0, n = 0;
  function r(s, o, a, l) {
    i = s, e = a, t = -3 * s + 3 * o - 2 * a - l, n = 2 * s - 2 * o + a + l;
  }
  return {
    initCatmullRom: function(s, o, a, l, c) {
      r(o, a, c * (a - s), c * (l - o));
    },
    initNonuniformCatmullRom: function(s, o, a, l, c, h, f) {
      let u = (o - s) / c - (a - s) / (c + h) + (a - o) / h, m = (a - o) / h - (l - o) / (h + f) + (l - a) / f;
      u *= h, m *= h, r(o, a, u, m);
    },
    calc: function(s) {
      const o = s * s, a = o * s;
      return i + e * s + t * o + n * a;
    }
  };
}
const Ki = /* @__PURE__ */ new U(), Gr = /* @__PURE__ */ new ps(), Vr = /* @__PURE__ */ new ps(), kr = /* @__PURE__ */ new ps();
class _p extends Bt {
  constructor(e = [], t = !1, n = "centripetal", r = 0.5) {
    super(), this.isCatmullRomCurve3 = !0, this.type = "CatmullRomCurve3", this.points = e, this.closed = t, this.curveType = n, this.tension = r;
  }
  getPoint(e, t = new U()) {
    const n = t, r = this.points, s = r.length, o = (s - (this.closed ? 0 : 1)) * e;
    let a = Math.floor(o), l = o - a;
    this.closed ? a += a > 0 ? 0 : (Math.floor(Math.abs(a) / s) + 1) * s : l === 0 && a === s - 1 && (a = s - 2, l = 1);
    let c, h;
    this.closed || a > 0 ? c = r[(a - 1) % s] : (Ki.subVectors(r[0], r[1]).add(r[0]), c = Ki);
    const f = r[a % s], u = r[(a + 1) % s];
    if (this.closed || a + 2 < s ? h = r[(a + 2) % s] : (Ki.subVectors(r[s - 1], r[s - 2]).add(r[s - 1]), h = Ki), this.curveType === "centripetal" || this.curveType === "chordal") {
      const m = this.curveType === "chordal" ? 0.5 : 0.25;
      let g = Math.pow(c.distanceToSquared(f), m), x = Math.pow(f.distanceToSquared(u), m), p = Math.pow(u.distanceToSquared(h), m);
      x < 1e-4 && (x = 1), g < 1e-4 && (g = x), p < 1e-4 && (p = x), Gr.initNonuniformCatmullRom(c.x, f.x, u.x, h.x, g, x, p), Vr.initNonuniformCatmullRom(c.y, f.y, u.y, h.y, g, x, p), kr.initNonuniformCatmullRom(c.z, f.z, u.z, h.z, g, x, p);
    } else
      this.curveType === "catmullrom" && (Gr.initCatmullRom(c.x, f.x, u.x, h.x, this.tension), Vr.initCatmullRom(c.y, f.y, u.y, h.y, this.tension), kr.initCatmullRom(c.z, f.z, u.z, h.z, this.tension));
    return n.set(
      Gr.calc(l),
      Vr.calc(l),
      kr.calc(l)
    ), n;
  }
  copy(e) {
    super.copy(e), this.points = [];
    for (let t = 0, n = e.points.length; t < n; t++) {
      const r = e.points[t];
      this.points.push(r.clone());
    }
    return this.closed = e.closed, this.curveType = e.curveType, this.tension = e.tension, this;
  }
  toJSON() {
    const e = super.toJSON();
    e.points = [];
    for (let t = 0, n = this.points.length; t < n; t++) {
      const r = this.points[t];
      e.points.push(r.toArray());
    }
    return e.closed = this.closed, e.curveType = this.curveType, e.tension = this.tension, e;
  }
  fromJSON(e) {
    super.fromJSON(e), this.points = [];
    for (let t = 0, n = e.points.length; t < n; t++) {
      const r = e.points[t];
      this.points.push(new U().fromArray(r));
    }
    return this.closed = e.closed, this.curveType = e.curveType, this.tension = e.tension, this;
  }
}
function za(i, e, t, n, r) {
  const s = (n - e) * 0.5, o = (r - t) * 0.5, a = i * i, l = i * a;
  return (2 * t - 2 * n + s + o) * l + (-3 * t + 3 * n - 2 * s - o) * a + s * i + t;
}
function vp(i, e) {
  const t = 1 - i;
  return t * t * e;
}
function xp(i, e) {
  return 2 * (1 - i) * i * e;
}
function Mp(i, e) {
  return i * i * e;
}
function pi(i, e, t, n) {
  return vp(i, e) + xp(i, t) + Mp(i, n);
}
function Sp(i, e) {
  const t = 1 - i;
  return t * t * t * e;
}
function Ep(i, e) {
  const t = 1 - i;
  return 3 * t * t * i * e;
}
function yp(i, e) {
  return 3 * (1 - i) * i * i * e;
}
function Tp(i, e) {
  return i * i * i * e;
}
function mi(i, e, t, n, r) {
  return Sp(i, e) + Ep(i, t) + yp(i, n) + Tp(i, r);
}
class Lo extends Bt {
  constructor(e = new oe(), t = new oe(), n = new oe(), r = new oe()) {
    super(), this.isCubicBezierCurve = !0, this.type = "CubicBezierCurve", this.v0 = e, this.v1 = t, this.v2 = n, this.v3 = r;
  }
  getPoint(e, t = new oe()) {
    const n = t, r = this.v0, s = this.v1, o = this.v2, a = this.v3;
    return n.set(
      mi(e, r.x, s.x, o.x, a.x),
      mi(e, r.y, s.y, o.y, a.y)
    ), n;
  }
  copy(e) {
    return super.copy(e), this.v0.copy(e.v0), this.v1.copy(e.v1), this.v2.copy(e.v2), this.v3.copy(e.v3), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v0 = this.v0.toArray(), e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e.v3 = this.v3.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v0.fromArray(e.v0), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this.v3.fromArray(e.v3), this;
  }
}
class Ap extends Bt {
  constructor(e = new U(), t = new U(), n = new U(), r = new U()) {
    super(), this.isCubicBezierCurve3 = !0, this.type = "CubicBezierCurve3", this.v0 = e, this.v1 = t, this.v2 = n, this.v3 = r;
  }
  getPoint(e, t = new U()) {
    const n = t, r = this.v0, s = this.v1, o = this.v2, a = this.v3;
    return n.set(
      mi(e, r.x, s.x, o.x, a.x),
      mi(e, r.y, s.y, o.y, a.y),
      mi(e, r.z, s.z, o.z, a.z)
    ), n;
  }
  copy(e) {
    return super.copy(e), this.v0.copy(e.v0), this.v1.copy(e.v1), this.v2.copy(e.v2), this.v3.copy(e.v3), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v0 = this.v0.toArray(), e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e.v3 = this.v3.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v0.fromArray(e.v0), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this.v3.fromArray(e.v3), this;
  }
}
class ms extends Bt {
  constructor(e = new oe(), t = new oe()) {
    super(), this.isLineCurve = !0, this.type = "LineCurve", this.v1 = e, this.v2 = t;
  }
  getPoint(e, t = new oe()) {
    const n = t;
    return e === 1 ? n.copy(this.v2) : (n.copy(this.v2).sub(this.v1), n.multiplyScalar(e).add(this.v1)), n;
  }
  // Line curve is linear, so we can overwrite default getPointAt
  getPointAt(e, t) {
    return this.getPoint(e, t);
  }
  getTangent(e, t = new oe()) {
    return t.subVectors(this.v2, this.v1).normalize();
  }
  getTangentAt(e, t) {
    return this.getTangent(e, t);
  }
  copy(e) {
    return super.copy(e), this.v1.copy(e.v1), this.v2.copy(e.v2), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this;
  }
}
class bp extends Bt {
  constructor(e = new U(), t = new U()) {
    super(), this.isLineCurve3 = !0, this.type = "LineCurve3", this.v1 = e, this.v2 = t;
  }
  getPoint(e, t = new U()) {
    const n = t;
    return e === 1 ? n.copy(this.v2) : (n.copy(this.v2).sub(this.v1), n.multiplyScalar(e).add(this.v1)), n;
  }
  // Line curve is linear, so we can overwrite default getPointAt
  getPointAt(e, t) {
    return this.getPoint(e, t);
  }
  getTangent(e, t = new U()) {
    return t.subVectors(this.v2, this.v1).normalize();
  }
  getTangentAt(e, t) {
    return this.getTangent(e, t);
  }
  copy(e) {
    return super.copy(e), this.v1.copy(e.v1), this.v2.copy(e.v2), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this;
  }
}
class Uo extends Bt {
  constructor(e = new oe(), t = new oe(), n = new oe()) {
    super(), this.isQuadraticBezierCurve = !0, this.type = "QuadraticBezierCurve", this.v0 = e, this.v1 = t, this.v2 = n;
  }
  getPoint(e, t = new oe()) {
    const n = t, r = this.v0, s = this.v1, o = this.v2;
    return n.set(
      pi(e, r.x, s.x, o.x),
      pi(e, r.y, s.y, o.y)
    ), n;
  }
  copy(e) {
    return super.copy(e), this.v0.copy(e.v0), this.v1.copy(e.v1), this.v2.copy(e.v2), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v0 = this.v0.toArray(), e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v0.fromArray(e.v0), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this;
  }
}
class wp extends Bt {
  constructor(e = new U(), t = new U(), n = new U()) {
    super(), this.isQuadraticBezierCurve3 = !0, this.type = "QuadraticBezierCurve3", this.v0 = e, this.v1 = t, this.v2 = n;
  }
  getPoint(e, t = new U()) {
    const n = t, r = this.v0, s = this.v1, o = this.v2;
    return n.set(
      pi(e, r.x, s.x, o.x),
      pi(e, r.y, s.y, o.y),
      pi(e, r.z, s.z, o.z)
    ), n;
  }
  copy(e) {
    return super.copy(e), this.v0.copy(e.v0), this.v1.copy(e.v1), this.v2.copy(e.v2), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.v0 = this.v0.toArray(), e.v1 = this.v1.toArray(), e.v2 = this.v2.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.v0.fromArray(e.v0), this.v1.fromArray(e.v1), this.v2.fromArray(e.v2), this;
  }
}
class Do extends Bt {
  constructor(e = []) {
    super(), this.isSplineCurve = !0, this.type = "SplineCurve", this.points = e;
  }
  getPoint(e, t = new oe()) {
    const n = t, r = this.points, s = (r.length - 1) * e, o = Math.floor(s), a = s - o, l = r[o === 0 ? o : o - 1], c = r[o], h = r[o > r.length - 2 ? r.length - 1 : o + 1], f = r[o > r.length - 3 ? r.length - 1 : o + 2];
    return n.set(
      za(a, l.x, c.x, h.x, f.x),
      za(a, l.y, c.y, h.y, f.y)
    ), n;
  }
  copy(e) {
    super.copy(e), this.points = [];
    for (let t = 0, n = e.points.length; t < n; t++) {
      const r = e.points[t];
      this.points.push(r.clone());
    }
    return this;
  }
  toJSON() {
    const e = super.toJSON();
    e.points = [];
    for (let t = 0, n = this.points.length; t < n; t++) {
      const r = this.points[t];
      e.points.push(r.toArray());
    }
    return e;
  }
  fromJSON(e) {
    super.fromJSON(e), this.points = [];
    for (let t = 0, n = e.points.length; t < n; t++) {
      const r = e.points[t];
      this.points.push(new oe().fromArray(r));
    }
    return this;
  }
}
var Io = /* @__PURE__ */ Object.freeze({
  __proto__: null,
  ArcCurve: gp,
  CatmullRomCurve3: _p,
  CubicBezierCurve: Lo,
  CubicBezierCurve3: Ap,
  EllipseCurve: ds,
  LineCurve: ms,
  LineCurve3: bp,
  QuadraticBezierCurve: Uo,
  QuadraticBezierCurve3: wp,
  SplineCurve: Do
});
class Rp extends Bt {
  constructor() {
    super(), this.type = "CurvePath", this.curves = [], this.autoClose = !1;
  }
  add(e) {
    this.curves.push(e);
  }
  closePath() {
    const e = this.curves[0].getPoint(0), t = this.curves[this.curves.length - 1].getPoint(1);
    e.equals(t) || this.curves.push(new ms(t, e));
  }
  // To get accurate point with reference to
  // entire path distance at time t,
  // following has to be done:
  // 1. Length of each sub path have to be known
  // 2. Locate and identify type of curve
  // 3. Get t for the curve
  // 4. Return curve.getPointAt(t')
  getPoint(e, t) {
    const n = e * this.getLength(), r = this.getCurveLengths();
    let s = 0;
    for (; s < r.length; ) {
      if (r[s] >= n) {
        const o = r[s] - n, a = this.curves[s], l = a.getLength(), c = l === 0 ? 0 : 1 - o / l;
        return a.getPointAt(c, t);
      }
      s++;
    }
    return null;
  }
  // We cannot use the default THREE.Curve getPoint() with getLength() because in
  // THREE.Curve, getLength() depends on getPoint() but in THREE.CurvePath
  // getPoint() depends on getLength
  getLength() {
    const e = this.getCurveLengths();
    return e[e.length - 1];
  }
  // cacheLengths must be recalculated.
  updateArcLengths() {
    this.needsUpdate = !0, this.cacheLengths = null, this.getCurveLengths();
  }
  // Compute lengths and cache them
  // We cannot overwrite getLengths() because UtoT mapping uses it.
  getCurveLengths() {
    if (this.cacheLengths && this.cacheLengths.length === this.curves.length)
      return this.cacheLengths;
    const e = [];
    let t = 0;
    for (let n = 0, r = this.curves.length; n < r; n++)
      t += this.curves[n].getLength(), e.push(t);
    return this.cacheLengths = e, e;
  }
  getSpacedPoints(e = 40) {
    const t = [];
    for (let n = 0; n <= e; n++)
      t.push(this.getPoint(n / e));
    return this.autoClose && t.push(t[0]), t;
  }
  getPoints(e = 12) {
    const t = [];
    let n;
    for (let r = 0, s = this.curves; r < s.length; r++) {
      const o = s[r], a = o.isEllipseCurve ? e * 2 : o.isLineCurve || o.isLineCurve3 ? 1 : o.isSplineCurve ? e * o.points.length : e, l = o.getPoints(a);
      for (let c = 0; c < l.length; c++) {
        const h = l[c];
        n && n.equals(h) || (t.push(h), n = h);
      }
    }
    return this.autoClose && t.length > 1 && !t[t.length - 1].equals(t[0]) && t.push(t[0]), t;
  }
  copy(e) {
    super.copy(e), this.curves = [];
    for (let t = 0, n = e.curves.length; t < n; t++) {
      const r = e.curves[t];
      this.curves.push(r.clone());
    }
    return this.autoClose = e.autoClose, this;
  }
  toJSON() {
    const e = super.toJSON();
    e.autoClose = this.autoClose, e.curves = [];
    for (let t = 0, n = this.curves.length; t < n; t++) {
      const r = this.curves[t];
      e.curves.push(r.toJSON());
    }
    return e;
  }
  fromJSON(e) {
    super.fromJSON(e), this.autoClose = e.autoClose, this.curves = [];
    for (let t = 0, n = e.curves.length; t < n; t++) {
      const r = e.curves[t];
      this.curves.push(new Io[r.type]().fromJSON(r));
    }
    return this;
  }
}
class Ha extends Rp {
  constructor(e) {
    super(), this.type = "Path", this.currentPoint = new oe(), e && this.setFromPoints(e);
  }
  setFromPoints(e) {
    this.moveTo(e[0].x, e[0].y);
    for (let t = 1, n = e.length; t < n; t++)
      this.lineTo(e[t].x, e[t].y);
    return this;
  }
  moveTo(e, t) {
    return this.currentPoint.set(e, t), this;
  }
  lineTo(e, t) {
    const n = new ms(this.currentPoint.clone(), new oe(e, t));
    return this.curves.push(n), this.currentPoint.set(e, t), this;
  }
  quadraticCurveTo(e, t, n, r) {
    const s = new Uo(
      this.currentPoint.clone(),
      new oe(e, t),
      new oe(n, r)
    );
    return this.curves.push(s), this.currentPoint.set(n, r), this;
  }
  bezierCurveTo(e, t, n, r, s, o) {
    const a = new Lo(
      this.currentPoint.clone(),
      new oe(e, t),
      new oe(n, r),
      new oe(s, o)
    );
    return this.curves.push(a), this.currentPoint.set(s, o), this;
  }
  splineThru(e) {
    const t = [this.currentPoint.clone()].concat(e), n = new Do(t);
    return this.curves.push(n), this.currentPoint.copy(e[e.length - 1]), this;
  }
  arc(e, t, n, r, s, o) {
    const a = this.currentPoint.x, l = this.currentPoint.y;
    return this.absarc(
      e + a,
      t + l,
      n,
      r,
      s,
      o
    ), this;
  }
  absarc(e, t, n, r, s, o) {
    return this.absellipse(e, t, n, n, r, s, o), this;
  }
  ellipse(e, t, n, r, s, o, a, l) {
    const c = this.currentPoint.x, h = this.currentPoint.y;
    return this.absellipse(e + c, t + h, n, r, s, o, a, l), this;
  }
  absellipse(e, t, n, r, s, o, a, l) {
    const c = new ds(e, t, n, r, s, o, a, l);
    if (this.curves.length > 0) {
      const f = c.getPoint(0);
      f.equals(this.currentPoint) || this.lineTo(f.x, f.y);
    }
    this.curves.push(c);
    const h = c.getPoint(1);
    return this.currentPoint.copy(h), this;
  }
  copy(e) {
    return super.copy(e), this.currentPoint.copy(e.currentPoint), this;
  }
  toJSON() {
    const e = super.toJSON();
    return e.currentPoint = this.currentPoint.toArray(), e;
  }
  fromJSON(e) {
    return super.fromJSON(e), this.currentPoint.fromArray(e.currentPoint), this;
  }
}
class No extends Ha {
  constructor(e) {
    super(e), this.uuid = Cn(), this.type = "Shape", this.holes = [];
  }
  getPointsHoles(e) {
    const t = [];
    for (let n = 0, r = this.holes.length; n < r; n++)
      t[n] = this.holes[n].getPoints(e);
    return t;
  }
  // get points of shape and holes (keypoints based on segments parameter)
  extractPoints(e) {
    return {
      shape: this.getPoints(e),
      holes: this.getPointsHoles(e)
    };
  }
  copy(e) {
    super.copy(e), this.holes = [];
    for (let t = 0, n = e.holes.length; t < n; t++) {
      const r = e.holes[t];
      this.holes.push(r.clone());
    }
    return this;
  }
  toJSON() {
    const e = super.toJSON();
    e.uuid = this.uuid, e.holes = [];
    for (let t = 0, n = this.holes.length; t < n; t++) {
      const r = this.holes[t];
      e.holes.push(r.toJSON());
    }
    return e;
  }
  fromJSON(e) {
    super.fromJSON(e), this.uuid = e.uuid, this.holes = [];
    for (let t = 0, n = e.holes.length; t < n; t++) {
      const r = e.holes[t];
      this.holes.push(new Ha().fromJSON(r));
    }
    return this;
  }
}
const Cp = {
  triangulate: function(i, e, t = 2) {
    const n = e && e.length, r = n ? e[0] * t : i.length;
    let s = Oo(i, 0, r, t, !0);
    const o = [];
    if (!s || s.next === s.prev)
      return o;
    let a, l, c, h, f, u, m;
    if (n && (s = Ip(i, e, s, t)), i.length > 80 * t) {
      a = c = i[0], l = h = i[1];
      for (let g = t; g < r; g += t)
        f = i[g], u = i[g + 1], f < a && (a = f), u < l && (l = u), f > c && (c = f), u > h && (h = u);
      m = Math.max(c - a, h - l), m = m !== 0 ? 32767 / m : 0;
    }
    return Si(s, o, t, a, l, m, 0), o;
  }
};
function Oo(i, e, t, n, r) {
  let s, o;
  if (r === Xp(i, e, t, n) > 0)
    for (s = e; s < t; s += n)
      o = Ga(s, i[s], i[s + 1], o);
  else
    for (s = t - n; s >= e; s -= n)
      o = Ga(s, i[s], i[s + 1], o);
  return o && lr(o, o.next) && (yi(o), o = o.next), o;
}
function wn(i, e) {
  if (!i)
    return i;
  e || (e = i);
  let t = i, n;
  do
    if (n = !1, !t.steiner && (lr(t, t.next) || Ke(t.prev, t, t.next) === 0)) {
      if (yi(t), t = e = t.prev, t === t.next)
        break;
      n = !0;
    } else
      t = t.next;
  while (n || t !== e);
  return e;
}
function Si(i, e, t, n, r, s, o) {
  if (!i)
    return;
  !o && s && zp(i, n, r, s);
  let a = i, l, c;
  for (; i.prev !== i.next; ) {
    if (l = i.prev, c = i.next, s ? Lp(i, n, r, s) : Pp(i)) {
      e.push(l.i / t | 0), e.push(i.i / t | 0), e.push(c.i / t | 0), yi(i), i = c.next, a = c.next;
      continue;
    }
    if (i = c, i === a) {
      o ? o === 1 ? (i = Up(wn(i), e, t), Si(i, e, t, n, r, s, 2)) : o === 2 && Dp(i, e, t, n, r, s) : Si(wn(i), e, t, n, r, s, 1);
      break;
    }
  }
}
function Pp(i) {
  const e = i.prev, t = i, n = i.next;
  if (Ke(e, t, n) >= 0)
    return !1;
  const r = e.x, s = t.x, o = n.x, a = e.y, l = t.y, c = n.y, h = r < s ? r < o ? r : o : s < o ? s : o, f = a < l ? a < c ? a : c : l < c ? l : c, u = r > s ? r > o ? r : o : s > o ? s : o, m = a > l ? a > c ? a : c : l > c ? l : c;
  let g = n.next;
  for (; g !== e; ) {
    if (g.x >= h && g.x <= u && g.y >= f && g.y <= m && Jn(r, a, s, l, o, c, g.x, g.y) && Ke(g.prev, g, g.next) >= 0)
      return !1;
    g = g.next;
  }
  return !0;
}
function Lp(i, e, t, n) {
  const r = i.prev, s = i, o = i.next;
  if (Ke(r, s, o) >= 0)
    return !1;
  const a = r.x, l = s.x, c = o.x, h = r.y, f = s.y, u = o.y, m = a < l ? a < c ? a : c : l < c ? l : c, g = h < f ? h < u ? h : u : f < u ? f : u, x = a > l ? a > c ? a : c : l > c ? l : c, p = h > f ? h > u ? h : u : f > u ? f : u, d = es(m, g, e, t, n), A = es(x, p, e, t, n);
  let _ = i.prevZ, T = i.nextZ;
  for (; _ && _.z >= d && T && T.z <= A; ) {
    if (_.x >= m && _.x <= x && _.y >= g && _.y <= p && _ !== r && _ !== o && Jn(a, h, l, f, c, u, _.x, _.y) && Ke(_.prev, _, _.next) >= 0 || (_ = _.prevZ, T.x >= m && T.x <= x && T.y >= g && T.y <= p && T !== r && T !== o && Jn(a, h, l, f, c, u, T.x, T.y) && Ke(T.prev, T, T.next) >= 0))
      return !1;
    T = T.nextZ;
  }
  for (; _ && _.z >= d; ) {
    if (_.x >= m && _.x <= x && _.y >= g && _.y <= p && _ !== r && _ !== o && Jn(a, h, l, f, c, u, _.x, _.y) && Ke(_.prev, _, _.next) >= 0)
      return !1;
    _ = _.prevZ;
  }
  for (; T && T.z <= A; ) {
    if (T.x >= m && T.x <= x && T.y >= g && T.y <= p && T !== r && T !== o && Jn(a, h, l, f, c, u, T.x, T.y) && Ke(T.prev, T, T.next) >= 0)
      return !1;
    T = T.nextZ;
  }
  return !0;
}
function Up(i, e, t) {
  let n = i;
  do {
    const r = n.prev, s = n.next.next;
    !lr(r, s) && Fo(r, n, n.next, s) && Ei(r, s) && Ei(s, r) && (e.push(r.i / t | 0), e.push(n.i / t | 0), e.push(s.i / t | 0), yi(n), yi(n.next), n = i = s), n = n.next;
  } while (n !== i);
  return wn(n);
}
function Dp(i, e, t, n, r, s) {
  let o = i;
  do {
    let a = o.next.next;
    for (; a !== o.prev; ) {
      if (o.i !== a.i && Vp(o, a)) {
        let l = Bo(o, a);
        o = wn(o, o.next), l = wn(l, l.next), Si(o, e, t, n, r, s, 0), Si(l, e, t, n, r, s, 0);
        return;
      }
      a = a.next;
    }
    o = o.next;
  } while (o !== i);
}
function Ip(i, e, t, n) {
  const r = [];
  let s, o, a, l, c;
  for (s = 0, o = e.length; s < o; s++)
    a = e[s] * n, l = s < o - 1 ? e[s + 1] * n : i.length, c = Oo(i, a, l, n, !1), c === c.next && (c.steiner = !0), r.push(Gp(c));
  for (r.sort(Np), s = 0; s < r.length; s++)
    t = Op(r[s], t);
  return t;
}
function Np(i, e) {
  return i.x - e.x;
}
function Op(i, e) {
  const t = Fp(i, e);
  if (!t)
    return e;
  const n = Bo(t, i);
  return wn(n, n.next), wn(t, t.next);
}
function Fp(i, e) {
  let t = e, n = -1 / 0, r;
  const s = i.x, o = i.y;
  do {
    if (o <= t.y && o >= t.next.y && t.next.y !== t.y) {
      const u = t.x + (o - t.y) * (t.next.x - t.x) / (t.next.y - t.y);
      if (u <= s && u > n && (n = u, r = t.x < t.next.x ? t : t.next, u === s))
        return r;
    }
    t = t.next;
  } while (t !== e);
  if (!r)
    return null;
  const a = r, l = r.x, c = r.y;
  let h = 1 / 0, f;
  t = r;
  do
    s >= t.x && t.x >= l && s !== t.x && Jn(o < c ? s : n, o, l, c, o < c ? n : s, o, t.x, t.y) && (f = Math.abs(o - t.y) / (s - t.x), Ei(t, i) && (f < h || f === h && (t.x > r.x || t.x === r.x && Bp(r, t))) && (r = t, h = f)), t = t.next;
  while (t !== a);
  return r;
}
function Bp(i, e) {
  return Ke(i.prev, i, e.prev) < 0 && Ke(e.next, i, i.next) < 0;
}
function zp(i, e, t, n) {
  let r = i;
  do
    r.z === 0 && (r.z = es(r.x, r.y, e, t, n)), r.prevZ = r.prev, r.nextZ = r.next, r = r.next;
  while (r !== i);
  r.prevZ.nextZ = null, r.prevZ = null, Hp(r);
}
function Hp(i) {
  let e, t, n, r, s, o, a, l, c = 1;
  do {
    for (t = i, i = null, s = null, o = 0; t; ) {
      for (o++, n = t, a = 0, e = 0; e < c && (a++, n = n.nextZ, !!n); e++)
        ;
      for (l = c; a > 0 || l > 0 && n; )
        a !== 0 && (l === 0 || !n || t.z <= n.z) ? (r = t, t = t.nextZ, a--) : (r = n, n = n.nextZ, l--), s ? s.nextZ = r : i = r, r.prevZ = s, s = r;
      t = n;
    }
    s.nextZ = null, c *= 2;
  } while (o > 1);
  return i;
}
function es(i, e, t, n, r) {
  return i = (i - t) * r | 0, e = (e - n) * r | 0, i = (i | i << 8) & 16711935, i = (i | i << 4) & 252645135, i = (i | i << 2) & 858993459, i = (i | i << 1) & 1431655765, e = (e | e << 8) & 16711935, e = (e | e << 4) & 252645135, e = (e | e << 2) & 858993459, e = (e | e << 1) & 1431655765, i | e << 1;
}
function Gp(i) {
  let e = i, t = i;
  do
    (e.x < t.x || e.x === t.x && e.y < t.y) && (t = e), e = e.next;
  while (e !== i);
  return t;
}
function Jn(i, e, t, n, r, s, o, a) {
  return (r - o) * (e - a) >= (i - o) * (s - a) && (i - o) * (n - a) >= (t - o) * (e - a) && (t - o) * (s - a) >= (r - o) * (n - a);
}
function Vp(i, e) {
  return i.next.i !== e.i && i.prev.i !== e.i && !kp(i, e) && // dones't intersect other edges
  (Ei(i, e) && Ei(e, i) && Wp(i, e) && // locally visible
  (Ke(i.prev, i, e.prev) || Ke(i, e.prev, e)) || // does not create opposite-facing sectors
  lr(i, e) && Ke(i.prev, i, i.next) > 0 && Ke(e.prev, e, e.next) > 0);
}
function Ke(i, e, t) {
  return (e.y - i.y) * (t.x - e.x) - (e.x - i.x) * (t.y - e.y);
}
function lr(i, e) {
  return i.x === e.x && i.y === e.y;
}
function Fo(i, e, t, n) {
  const r = $i(Ke(i, e, t)), s = $i(Ke(i, e, n)), o = $i(Ke(t, n, i)), a = $i(Ke(t, n, e));
  return !!(r !== s && o !== a || r === 0 && Ji(i, t, e) || s === 0 && Ji(i, n, e) || o === 0 && Ji(t, i, n) || a === 0 && Ji(t, e, n));
}
function Ji(i, e, t) {
  return e.x <= Math.max(i.x, t.x) && e.x >= Math.min(i.x, t.x) && e.y <= Math.max(i.y, t.y) && e.y >= Math.min(i.y, t.y);
}
function $i(i) {
  return i > 0 ? 1 : i < 0 ? -1 : 0;
}
function kp(i, e) {
  let t = i;
  do {
    if (t.i !== i.i && t.next.i !== i.i && t.i !== e.i && t.next.i !== e.i && Fo(t, t.next, i, e))
      return !0;
    t = t.next;
  } while (t !== i);
  return !1;
}
function Ei(i, e) {
  return Ke(i.prev, i, i.next) < 0 ? Ke(i, e, i.next) >= 0 && Ke(i, i.prev, e) >= 0 : Ke(i, e, i.prev) < 0 || Ke(i, i.next, e) < 0;
}
function Wp(i, e) {
  let t = i, n = !1;
  const r = (i.x + e.x) / 2, s = (i.y + e.y) / 2;
  do
    t.y > s != t.next.y > s && t.next.y !== t.y && r < (t.next.x - t.x) * (s - t.y) / (t.next.y - t.y) + t.x && (n = !n), t = t.next;
  while (t !== i);
  return n;
}
function Bo(i, e) {
  const t = new ts(i.i, i.x, i.y), n = new ts(e.i, e.x, e.y), r = i.next, s = e.prev;
  return i.next = e, e.prev = i, t.next = r, r.prev = t, n.next = t, t.prev = n, s.next = n, n.prev = s, n;
}
function Ga(i, e, t, n) {
  const r = new ts(i, e, t);
  return n ? (r.next = n.next, r.prev = n, n.next.prev = r, n.next = r) : (r.prev = r, r.next = r), r;
}
function yi(i) {
  i.next.prev = i.prev, i.prev.next = i.next, i.prevZ && (i.prevZ.nextZ = i.nextZ), i.nextZ && (i.nextZ.prevZ = i.prevZ);
}
function ts(i, e, t) {
  this.i = i, this.x = e, this.y = t, this.prev = null, this.next = null, this.z = 0, this.prevZ = null, this.nextZ = null, this.steiner = !1;
}
function Xp(i, e, t, n) {
  let r = 0;
  for (let s = e, o = t - n; s < t; s += n)
    r += (i[o] - i[s]) * (i[s + 1] + i[o + 1]), o = s;
  return r;
}
class gi {
  // calculate area of the contour polygon
  static area(e) {
    const t = e.length;
    let n = 0;
    for (let r = t - 1, s = 0; s < t; r = s++)
      n += e[r].x * e[s].y - e[s].x * e[r].y;
    return n * 0.5;
  }
  static isClockWise(e) {
    return gi.area(e) < 0;
  }
  static triangulateShape(e, t) {
    const n = [], r = [], s = [];
    Va(e), ka(n, e);
    let o = e.length;
    t.forEach(Va);
    for (let l = 0; l < t.length; l++)
      r.push(o), o += t[l].length, ka(n, t[l]);
    const a = Cp.triangulate(n, r);
    for (let l = 0; l < a.length; l += 3)
      s.push(a.slice(l, l + 3));
    return s;
  }
}
function Va(i) {
  const e = i.length;
  e > 2 && i[e - 1].equals(i[0]) && i.pop();
}
function ka(i, e) {
  for (let t = 0; t < e.length; t++)
    i.push(e[t].x), i.push(e[t].y);
}
class gs extends cn {
  constructor(e = new No([new oe(0.5, 0.5), new oe(-0.5, 0.5), new oe(-0.5, -0.5), new oe(0.5, -0.5)]), t = {}) {
    super(), this.type = "ExtrudeGeometry", this.parameters = {
      shapes: e,
      options: t
    }, e = Array.isArray(e) ? e : [e];
    const n = this, r = [], s = [];
    for (let a = 0, l = e.length; a < l; a++) {
      const c = e[a];
      o(c);
    }
    this.setAttribute("position", new jt(r, 3)), this.setAttribute("uv", new jt(s, 2)), this.computeVertexNormals();
    function o(a) {
      const l = [], c = t.curveSegments !== void 0 ? t.curveSegments : 12, h = t.steps !== void 0 ? t.steps : 1, f = t.depth !== void 0 ? t.depth : 1;
      let u = t.bevelEnabled !== void 0 ? t.bevelEnabled : !0, m = t.bevelThickness !== void 0 ? t.bevelThickness : 0.2, g = t.bevelSize !== void 0 ? t.bevelSize : m - 0.1, x = t.bevelOffset !== void 0 ? t.bevelOffset : 0, p = t.bevelSegments !== void 0 ? t.bevelSegments : 3;
      const d = t.extrudePath, A = t.UVGenerator !== void 0 ? t.UVGenerator : Yp;
      let _, T = !1, C, L, w, V;
      d && (_ = d.getSpacedPoints(h), T = !0, u = !1, C = d.computeFrenetFrames(h, !1), L = new U(), w = new U(), V = new U()), u || (p = 0, m = 0, g = 0, x = 0);
      const M = a.extractPoints(c);
      let y = M.shape;
      const Y = M.holes;
      if (!gi.isClockWise(y)) {
        y = y.reverse();
        for (let b = 0, le = Y.length; b < le; b++) {
          const Z = Y[b];
          gi.isClockWise(Z) && (Y[b] = Z.reverse());
        }
      }
      const B = gi.triangulateShape(y, Y), H = y;
      for (let b = 0, le = Y.length; b < le; b++) {
        const Z = Y[b];
        y = y.concat(Z);
      }
      function G(b, le, Z) {
        return le || console.error("THREE.ExtrudeGeometry: vec does not exist"), b.clone().addScaledVector(le, Z);
      }
      const Q = y.length, X = B.length;
      function j(b, le, Z) {
        let re, $, ye;
        const xe = b.x - le.x, Me = b.y - le.y, Ce = Z.x - b.x, ze = Z.y - b.y, Ze = xe * xe + Me * Me, E = xe * ze - Me * Ce;
        if (Math.abs(E) > Number.EPSILON) {
          const v = Math.sqrt(Ze), O = Math.sqrt(Ce * Ce + ze * ze), ie = le.x - Me / v, te = le.y + xe / v, se = Z.x - ze / O, Ee = Z.y + Ce / O, ae = ((se - ie) * ze - (Ee - te) * Ce) / (xe * ze - Me * Ce);
          re = ie + xe * ae - b.x, $ = te + Me * ae - b.y;
          const z = re * re + $ * $;
          if (z <= 2)
            return new oe(re, $);
          ye = Math.sqrt(z / 2);
        } else {
          let v = !1;
          xe > Number.EPSILON ? Ce > Number.EPSILON && (v = !0) : xe < -Number.EPSILON ? Ce < -Number.EPSILON && (v = !0) : Math.sign(Me) === Math.sign(ze) && (v = !0), v ? (re = -Me, $ = xe, ye = Math.sqrt(Ze)) : (re = xe, $ = Me, ye = Math.sqrt(Ze / 2));
        }
        return new oe(re / ye, $ / ye);
      }
      const J = [];
      for (let b = 0, le = H.length, Z = le - 1, re = b + 1; b < le; b++, Z++, re++)
        Z === le && (Z = 0), re === le && (re = 0), J[b] = j(H[b], H[Z], H[re]);
      const ee = [];
      let I, q = J.concat();
      for (let b = 0, le = Y.length; b < le; b++) {
        const Z = Y[b];
        I = [];
        for (let re = 0, $ = Z.length, ye = $ - 1, xe = re + 1; re < $; re++, ye++, xe++)
          ye === $ && (ye = 0), xe === $ && (xe = 0), I[re] = j(Z[re], Z[ye], Z[xe]);
        ee.push(I), q = q.concat(I);
      }
      for (let b = 0; b < p; b++) {
        const le = b / p, Z = m * Math.cos(le * Math.PI / 2), re = g * Math.sin(le * Math.PI / 2) + x;
        for (let $ = 0, ye = H.length; $ < ye; $++) {
          const xe = G(H[$], J[$], re);
          be(xe.x, xe.y, -Z);
        }
        for (let $ = 0, ye = Y.length; $ < ye; $++) {
          const xe = Y[$];
          I = ee[$];
          for (let Me = 0, Ce = xe.length; Me < Ce; Me++) {
            const ze = G(xe[Me], I[Me], re);
            be(ze.x, ze.y, -Z);
          }
        }
      }
      const pe = g + x;
      for (let b = 0; b < Q; b++) {
        const le = u ? G(y[b], q[b], pe) : y[b];
        T ? (w.copy(C.normals[0]).multiplyScalar(le.x), L.copy(C.binormals[0]).multiplyScalar(le.y), V.copy(_[0]).add(w).add(L), be(V.x, V.y, V.z)) : be(le.x, le.y, 0);
      }
      for (let b = 1; b <= h; b++)
        for (let le = 0; le < Q; le++) {
          const Z = u ? G(y[le], q[le], pe) : y[le];
          T ? (w.copy(C.normals[b]).multiplyScalar(Z.x), L.copy(C.binormals[b]).multiplyScalar(Z.y), V.copy(_[b]).add(w).add(L), be(V.x, V.y, V.z)) : be(Z.x, Z.y, f / h * b);
        }
      for (let b = p - 1; b >= 0; b--) {
        const le = b / p, Z = m * Math.cos(le * Math.PI / 2), re = g * Math.sin(le * Math.PI / 2) + x;
        for (let $ = 0, ye = H.length; $ < ye; $++) {
          const xe = G(H[$], J[$], re);
          be(xe.x, xe.y, f + Z);
        }
        for (let $ = 0, ye = Y.length; $ < ye; $++) {
          const xe = Y[$];
          I = ee[$];
          for (let Me = 0, Ce = xe.length; Me < Ce; Me++) {
            const ze = G(xe[Me], I[Me], re);
            T ? be(ze.x, ze.y + _[h - 1].y, _[h - 1].x + Z) : be(ze.x, ze.y, f + Z);
          }
        }
      }
      me(), ve();
      function me() {
        const b = r.length / 3;
        if (u) {
          let le = 0, Z = Q * le;
          for (let re = 0; re < X; re++) {
            const $ = B[re];
            Te($[2] + Z, $[1] + Z, $[0] + Z);
          }
          le = h + p * 2, Z = Q * le;
          for (let re = 0; re < X; re++) {
            const $ = B[re];
            Te($[0] + Z, $[1] + Z, $[2] + Z);
          }
        } else {
          for (let le = 0; le < X; le++) {
            const Z = B[le];
            Te(Z[2], Z[1], Z[0]);
          }
          for (let le = 0; le < X; le++) {
            const Z = B[le];
            Te(Z[0] + Q * h, Z[1] + Q * h, Z[2] + Q * h);
          }
        }
        n.addGroup(b, r.length / 3 - b, 0);
      }
      function ve() {
        const b = r.length / 3;
        let le = 0;
        Ae(H, le), le += H.length;
        for (let Z = 0, re = Y.length; Z < re; Z++) {
          const $ = Y[Z];
          Ae($, le), le += $.length;
        }
        n.addGroup(b, r.length / 3 - b, 1);
      }
      function Ae(b, le) {
        let Z = b.length;
        for (; --Z >= 0; ) {
          const re = Z;
          let $ = Z - 1;
          $ < 0 && ($ = b.length - 1);
          for (let ye = 0, xe = h + p * 2; ye < xe; ye++) {
            const Me = Q * ye, Ce = Q * (ye + 1), ze = le + re + Me, Ze = le + $ + Me, E = le + $ + Ce, v = le + re + Ce;
            ke(ze, Ze, E, v);
          }
        }
      }
      function be(b, le, Z) {
        l.push(b), l.push(le), l.push(Z);
      }
      function Te(b, le, Z) {
        Ye(b), Ye(le), Ye(Z);
        const re = r.length / 3, $ = A.generateTopUV(n, r, re - 3, re - 2, re - 1);
        we($[0]), we($[1]), we($[2]);
      }
      function ke(b, le, Z, re) {
        Ye(b), Ye(le), Ye(re), Ye(le), Ye(Z), Ye(re);
        const $ = r.length / 3, ye = A.generateSideWallUV(n, r, $ - 6, $ - 3, $ - 2, $ - 1);
        we(ye[0]), we(ye[1]), we(ye[3]), we(ye[1]), we(ye[2]), we(ye[3]);
      }
      function Ye(b) {
        r.push(l[b * 3 + 0]), r.push(l[b * 3 + 1]), r.push(l[b * 3 + 2]);
      }
      function we(b) {
        s.push(b.x), s.push(b.y);
      }
    }
  }
  copy(e) {
    return super.copy(e), this.parameters = Object.assign({}, e.parameters), this;
  }
  toJSON() {
    const e = super.toJSON(), t = this.parameters.shapes, n = this.parameters.options;
    return qp(t, n, e);
  }
  static fromJSON(e, t) {
    const n = [];
    for (let s = 0, o = e.shapes.length; s < o; s++) {
      const a = t[e.shapes[s]];
      n.push(a);
    }
    const r = e.options.extrudePath;
    return r !== void 0 && (e.options.extrudePath = new Io[r.type]().fromJSON(r)), new gs(n, e.options);
  }
}
const Yp = {
  generateTopUV: function(i, e, t, n, r) {
    const s = e[t * 3], o = e[t * 3 + 1], a = e[n * 3], l = e[n * 3 + 1], c = e[r * 3], h = e[r * 3 + 1];
    return [
      new oe(s, o),
      new oe(a, l),
      new oe(c, h)
    ];
  },
  generateSideWallUV: function(i, e, t, n, r, s) {
    const o = e[t * 3], a = e[t * 3 + 1], l = e[t * 3 + 2], c = e[n * 3], h = e[n * 3 + 1], f = e[n * 3 + 2], u = e[r * 3], m = e[r * 3 + 1], g = e[r * 3 + 2], x = e[s * 3], p = e[s * 3 + 1], d = e[s * 3 + 2];
    return Math.abs(a - h) < Math.abs(o - c) ? [
      new oe(o, 1 - l),
      new oe(c, 1 - f),
      new oe(u, 1 - g),
      new oe(x, 1 - d)
    ] : [
      new oe(a, 1 - l),
      new oe(h, 1 - f),
      new oe(m, 1 - g),
      new oe(p, 1 - d)
    ];
  }
};
function qp(i, e, t) {
  if (t.shapes = [], Array.isArray(i))
    for (let n = 0, r = i.length; n < r; n++) {
      const s = i[n];
      t.shapes.push(s.uuid);
    }
  else
    t.shapes.push(i.uuid);
  return t.options = Object.assign({}, e), e.extrudePath !== void 0 && (t.options.extrudePath = e.extrudePath.toJSON()), t;
}
class zo extends Ai {
  constructor(e) {
    super(), this.isMeshStandardMaterial = !0, this.defines = { STANDARD: "" }, this.type = "MeshStandardMaterial", this.color = new We(16777215), this.roughness = 1, this.metalness = 0, this.map = null, this.lightMap = null, this.lightMapIntensity = 1, this.aoMap = null, this.aoMapIntensity = 1, this.emissive = new We(0), this.emissiveIntensity = 1, this.emissiveMap = null, this.bumpMap = null, this.bumpScale = 1, this.normalMap = null, this.normalMapType = ho, this.normalScale = new oe(1, 1), this.displacementMap = null, this.displacementScale = 1, this.displacementBias = 0, this.roughnessMap = null, this.metalnessMap = null, this.alphaMap = null, this.envMap = null, this.envMapIntensity = 1, this.wireframe = !1, this.wireframeLinewidth = 1, this.wireframeLinecap = "round", this.wireframeLinejoin = "round", this.flatShading = !1, this.fog = !0, this.setValues(e);
  }
  copy(e) {
    return super.copy(e), this.defines = { STANDARD: "" }, this.color.copy(e.color), this.roughness = e.roughness, this.metalness = e.metalness, this.map = e.map, this.lightMap = e.lightMap, this.lightMapIntensity = e.lightMapIntensity, this.aoMap = e.aoMap, this.aoMapIntensity = e.aoMapIntensity, this.emissive.copy(e.emissive), this.emissiveMap = e.emissiveMap, this.emissiveIntensity = e.emissiveIntensity, this.bumpMap = e.bumpMap, this.bumpScale = e.bumpScale, this.normalMap = e.normalMap, this.normalMapType = e.normalMapType, this.normalScale.copy(e.normalScale), this.displacementMap = e.displacementMap, this.displacementScale = e.displacementScale, this.displacementBias = e.displacementBias, this.roughnessMap = e.roughnessMap, this.metalnessMap = e.metalnessMap, this.alphaMap = e.alphaMap, this.envMap = e.envMap, this.envMapIntensity = e.envMapIntensity, this.wireframe = e.wireframe, this.wireframeLinewidth = e.wireframeLinewidth, this.wireframeLinecap = e.wireframeLinecap, this.wireframeLinejoin = e.wireframeLinejoin, this.flatShading = e.flatShading, this.fog = e.fog, this;
  }
}
class ns extends zo {
  constructor(e) {
    super(), this.isMeshPhysicalMaterial = !0, this.defines = {
      STANDARD: "",
      PHYSICAL: ""
    }, this.type = "MeshPhysicalMaterial", this.anisotropyRotation = 0, this.anisotropyMap = null, this.clearcoatMap = null, this.clearcoatRoughness = 0, this.clearcoatRoughnessMap = null, this.clearcoatNormalScale = new oe(1, 1), this.clearcoatNormalMap = null, this.ior = 1.5, Object.defineProperty(this, "reflectivity", {
      get: function() {
        return it(2.5 * (this.ior - 1) / (this.ior + 1), 0, 1);
      },
      set: function(t) {
        this.ior = (1 + 0.4 * t) / (1 - 0.4 * t);
      }
    }), this.iridescenceMap = null, this.iridescenceIOR = 1.3, this.iridescenceThicknessRange = [100, 400], this.iridescenceThicknessMap = null, this.sheenColor = new We(0), this.sheenColorMap = null, this.sheenRoughness = 1, this.sheenRoughnessMap = null, this.transmissionMap = null, this.thickness = 0, this.thicknessMap = null, this.attenuationDistance = 1 / 0, this.attenuationColor = new We(1, 1, 1), this.specularIntensity = 1, this.specularIntensityMap = null, this.specularColor = new We(1, 1, 1), this.specularColorMap = null, this._anisotropy = 0, this._clearcoat = 0, this._iridescence = 0, this._sheen = 0, this._transmission = 0, this.setValues(e);
  }
  get anisotropy() {
    return this._anisotropy;
  }
  set anisotropy(e) {
    this._anisotropy > 0 != e > 0 && this.version++, this._anisotropy = e;
  }
  get clearcoat() {
    return this._clearcoat;
  }
  set clearcoat(e) {
    this._clearcoat > 0 != e > 0 && this.version++, this._clearcoat = e;
  }
  get iridescence() {
    return this._iridescence;
  }
  set iridescence(e) {
    this._iridescence > 0 != e > 0 && this.version++, this._iridescence = e;
  }
  get sheen() {
    return this._sheen;
  }
  set sheen(e) {
    this._sheen > 0 != e > 0 && this.version++, this._sheen = e;
  }
  get transmission() {
    return this._transmission;
  }
  set transmission(e) {
    this._transmission > 0 != e > 0 && this.version++, this._transmission = e;
  }
  copy(e) {
    return super.copy(e), this.defines = {
      STANDARD: "",
      PHYSICAL: ""
    }, this.anisotropy = e.anisotropy, this.anisotropyRotation = e.anisotropyRotation, this.anisotropyMap = e.anisotropyMap, this.clearcoat = e.clearcoat, this.clearcoatMap = e.clearcoatMap, this.clearcoatRoughness = e.clearcoatRoughness, this.clearcoatRoughnessMap = e.clearcoatRoughnessMap, this.clearcoatNormalMap = e.clearcoatNormalMap, this.clearcoatNormalScale.copy(e.clearcoatNormalScale), this.ior = e.ior, this.iridescence = e.iridescence, this.iridescenceMap = e.iridescenceMap, this.iridescenceIOR = e.iridescenceIOR, this.iridescenceThicknessRange = [...e.iridescenceThicknessRange], this.iridescenceThicknessMap = e.iridescenceThicknessMap, this.sheen = e.sheen, this.sheenColor.copy(e.sheenColor), this.sheenColorMap = e.sheenColorMap, this.sheenRoughness = e.sheenRoughness, this.sheenRoughnessMap = e.sheenRoughnessMap, this.transmission = e.transmission, this.transmissionMap = e.transmissionMap, this.thickness = e.thickness, this.thicknessMap = e.thicknessMap, this.attenuationDistance = e.attenuationDistance, this.attenuationColor.copy(e.attenuationColor), this.specularIntensity = e.specularIntensity, this.specularIntensityMap = e.specularIntensityMap, this.specularColor.copy(e.specularColor), this.specularColorMap = e.specularColorMap, this;
  }
}
class Ho extends ht {
  constructor(e, t = 1) {
    super(), this.isLight = !0, this.type = "Light", this.color = new We(e), this.intensity = t;
  }
  dispose() {
  }
  copy(e, t) {
    return super.copy(e, t), this.color.copy(e.color), this.intensity = e.intensity, this;
  }
  toJSON(e) {
    const t = super.toJSON(e);
    return t.object.color = this.color.getHex(), t.object.intensity = this.intensity, this.groundColor !== void 0 && (t.object.groundColor = this.groundColor.getHex()), this.distance !== void 0 && (t.object.distance = this.distance), this.angle !== void 0 && (t.object.angle = this.angle), this.decay !== void 0 && (t.object.decay = this.decay), this.penumbra !== void 0 && (t.object.penumbra = this.penumbra), this.shadow !== void 0 && (t.object.shadow = this.shadow.toJSON()), t;
  }
}
const Wr = /* @__PURE__ */ new nt(), Wa = /* @__PURE__ */ new U(), Xa = /* @__PURE__ */ new U();
class jp {
  constructor(e) {
    this.camera = e, this.bias = 0, this.normalBias = 0, this.radius = 1, this.blurSamples = 8, this.mapSize = new oe(512, 512), this.map = null, this.mapPass = null, this.matrix = new nt(), this.autoUpdate = !0, this.needsUpdate = !1, this._frustum = new us(), this._frameExtents = new oe(1, 1), this._viewportCount = 1, this._viewports = [
      new ot(0, 0, 1, 1)
    ];
  }
  getViewportCount() {
    return this._viewportCount;
  }
  getFrustum() {
    return this._frustum;
  }
  updateMatrices(e) {
    const t = this.camera, n = this.matrix;
    Wa.setFromMatrixPosition(e.matrixWorld), t.position.copy(Wa), Xa.setFromMatrixPosition(e.target.matrixWorld), t.lookAt(Xa), t.updateMatrixWorld(), Wr.multiplyMatrices(t.projectionMatrix, t.matrixWorldInverse), this._frustum.setFromProjectionMatrix(Wr), n.set(
      0.5,
      0,
      0,
      0.5,
      0,
      0.5,
      0,
      0.5,
      0,
      0,
      0.5,
      0.5,
      0,
      0,
      0,
      1
    ), n.multiply(Wr);
  }
  getViewport(e) {
    return this._viewports[e];
  }
  getFrameExtents() {
    return this._frameExtents;
  }
  dispose() {
    this.map && this.map.dispose(), this.mapPass && this.mapPass.dispose();
  }
  copy(e) {
    return this.camera = e.camera.clone(), this.bias = e.bias, this.radius = e.radius, this.mapSize.copy(e.mapSize), this;
  }
  clone() {
    return new this.constructor().copy(this);
  }
  toJSON() {
    const e = {};
    return this.bias !== 0 && (e.bias = this.bias), this.normalBias !== 0 && (e.normalBias = this.normalBias), this.radius !== 1 && (e.radius = this.radius), (this.mapSize.x !== 512 || this.mapSize.y !== 512) && (e.mapSize = this.mapSize.toArray()), e.camera = this.camera.toJSON(!1).object, delete e.camera.matrix, e;
  }
}
class Zp extends jp {
  constructor() {
    super(new Ao(-5, 5, 5, -5, 0.5, 500)), this.isDirectionalLightShadow = !0;
  }
}
class Kp extends Ho {
  constructor(e, t) {
    super(e, t), this.isDirectionalLight = !0, this.type = "DirectionalLight", this.position.copy(ht.DEFAULT_UP), this.updateMatrix(), this.target = new ht(), this.shadow = new Zp();
  }
  dispose() {
    this.shadow.dispose();
  }
  copy(e) {
    return super.copy(e), this.target = e.target.clone(), this.shadow = e.shadow.clone(), this;
  }
}
class Jp extends Ho {
  constructor(e, t) {
    super(e, t), this.isAmbientLight = !0, this.type = "AmbientLight";
  }
}
class $p {
  constructor(e, t, n = 0, r = 1 / 0) {
    this.ray = new cs(e, t), this.near = n, this.far = r, this.camera = null, this.layers = new hs(), this.params = {
      Mesh: {},
      Line: { threshold: 1 },
      LOD: {},
      Points: { threshold: 1 },
      Sprite: {}
    };
  }
  set(e, t) {
    this.ray.set(e, t);
  }
  setFromCamera(e, t) {
    t.isPerspectiveCamera ? (this.ray.origin.setFromMatrixPosition(t.matrixWorld), this.ray.direction.set(e.x, e.y, 0.5).unproject(t).sub(this.ray.origin).normalize(), this.camera = t) : t.isOrthographicCamera ? (this.ray.origin.set(e.x, e.y, (t.near + t.far) / (t.near - t.far)).unproject(t), this.ray.direction.set(0, 0, -1).transformDirection(t.matrixWorld), this.camera = t) : console.error("THREE.Raycaster: Unsupported camera type: " + t.type);
  }
  intersectObject(e, t = !0, n = []) {
    return is(e, this, n, t), n.sort(Ya), n;
  }
  intersectObjects(e, t = !0, n = []) {
    for (let r = 0, s = e.length; r < s; r++)
      is(e[r], this, n, t);
    return n.sort(Ya), n;
  }
}
function Ya(i, e) {
  return i.distance - e.distance;
}
function is(i, e, t, n) {
  if (i.layers.test(e.layers) && i.raycast(e, t), n === !0) {
    const r = i.children;
    for (let s = 0, o = r.length; s < o; s++)
      is(r[s], e, t, !0);
  }
}
class qa {
  constructor(e = 1, t = 0, n = 0) {
    return this.radius = e, this.phi = t, this.theta = n, this;
  }
  set(e, t, n) {
    return this.radius = e, this.phi = t, this.theta = n, this;
  }
  copy(e) {
    return this.radius = e.radius, this.phi = e.phi, this.theta = e.theta, this;
  }
  // restrict phi to be between EPS and PI-EPS
  makeSafe() {
    return this.phi = Math.max(1e-6, Math.min(Math.PI - 1e-6, this.phi)), this;
  }
  setFromVector3(e) {
    return this.setFromCartesianCoords(e.x, e.y, e.z);
  }
  setFromCartesianCoords(e, t, n) {
    return this.radius = Math.sqrt(e * e + t * t + n * n), this.radius === 0 ? (this.theta = 0, this.phi = 0) : (this.theta = Math.atan2(e, n), this.phi = Math.acos(it(t / this.radius, -1, 1))), this;
  }
  clone() {
    return new this.constructor().copy(this);
  }
}
typeof __THREE_DEVTOOLS__ < "u" && __THREE_DEVTOOLS__.dispatchEvent(new CustomEvent("register", { detail: {
  revision: ss
} }));
typeof window < "u" && (window.__THREE__ ? console.warn("WARNING: Multiple instances of Three.js being imported.") : window.__THREE__ = ss);
const ja = { type: "change" }, Xr = { type: "start" }, Za = { type: "end" }, Qi = new cs(), Ka = new en(), Qp = Math.cos(70 * ac.DEG2RAD);
class em extends Rn {
  constructor(e, t) {
    super(), this.object = e, this.domElement = t, this.domElement.style.touchAction = "none", this.enabled = !0, this.target = new U(), this.minDistance = 0, this.maxDistance = 1 / 0, this.minZoom = 0, this.maxZoom = 1 / 0, this.minPolarAngle = 0, this.maxPolarAngle = Math.PI, this.minAzimuthAngle = -1 / 0, this.maxAzimuthAngle = 1 / 0, this.enableDamping = !1, this.dampingFactor = 0.05, this.enableZoom = !0, this.zoomSpeed = 1, this.enableRotate = !0, this.rotateSpeed = 1, this.enablePan = !0, this.panSpeed = 1, this.screenSpacePanning = !0, this.keyPanSpeed = 7, this.zoomToCursor = !1, this.autoRotate = !1, this.autoRotateSpeed = 2, this.keys = { LEFT: "ArrowLeft", UP: "ArrowUp", RIGHT: "ArrowRight", BOTTOM: "ArrowDown" }, this.mouseButtons = { LEFT: Pn.ROTATE, MIDDLE: Pn.DOLLY, RIGHT: Pn.PAN }, this.touches = { ONE: Ln.ROTATE, TWO: Ln.DOLLY_PAN }, this.target0 = this.target.clone(), this.position0 = this.object.position.clone(), this.zoom0 = this.object.zoom, this._domElementKeyEvents = null, this.getPolarAngle = function() {
      return a.phi;
    }, this.getAzimuthalAngle = function() {
      return a.theta;
    }, this.getDistance = function() {
      return this.object.position.distanceTo(this.target);
    }, this.listenToKeyEvents = function(R) {
      R.addEventListener("keydown", v), this._domElementKeyEvents = R;
    }, this.stopListenToKeyEvents = function() {
      this._domElementKeyEvents.removeEventListener("keydown", v), this._domElementKeyEvents = null;
    }, this.saveState = function() {
      n.target0.copy(n.target), n.position0.copy(n.object.position), n.zoom0 = n.object.zoom;
    }, this.reset = function() {
      n.target.copy(n.target0), n.object.position.copy(n.position0), n.object.zoom = n.zoom0, n.object.updateProjectionMatrix(), n.dispatchEvent(ja), n.update(), s = r.NONE;
    }, this.update = function() {
      const R = new U(), K = new An().setFromUnitVectors(e.up, new U(0, 1, 0)), _e = K.clone().invert(), he = new U(), ge = new An(), De = new U(), Ve = 2 * Math.PI;
      return function() {
        const fe = n.object.position;
        R.copy(fe).sub(n.target), R.applyQuaternion(K), a.setFromVector3(R), n.autoRotate && s === r.NONE && Y(M()), n.enableDamping ? (a.theta += l.theta * n.dampingFactor, a.phi += l.phi * n.dampingFactor) : (a.theta += l.theta, a.phi += l.phi);
        let F = n.minAzimuthAngle, ne = n.maxAzimuthAngle;
        isFinite(F) && isFinite(ne) && (F < -Math.PI ? F += Ve : F > Math.PI && (F -= Ve), ne < -Math.PI ? ne += Ve : ne > Math.PI && (ne -= Ve), F <= ne ? a.theta = Math.max(F, Math.min(ne, a.theta)) : a.theta = a.theta > (F + ne) / 2 ? Math.max(F, a.theta) : Math.min(ne, a.theta)), a.phi = Math.max(n.minPolarAngle, Math.min(n.maxPolarAngle, a.phi)), a.makeSafe(), n.enableDamping === !0 ? n.target.addScaledVector(h, n.dampingFactor) : n.target.add(h), n.zoomToCursor && L || n.object.isOrthographicCamera ? a.radius = J(a.radius) : a.radius = J(a.radius * c), R.setFromSpherical(a), R.applyQuaternion(_e), fe.copy(n.target).add(R), n.object.lookAt(n.target), n.enableDamping === !0 ? (l.theta *= 1 - n.dampingFactor, l.phi *= 1 - n.dampingFactor, h.multiplyScalar(1 - n.dampingFactor)) : (l.set(0, 0, 0), h.set(0, 0, 0));
        let de = !1;
        if (n.zoomToCursor && L) {
          let Fe = null;
          if (n.object.isPerspectiveCamera) {
            const Xe = R.length();
            Fe = J(Xe * c);
            const qe = Xe - Fe;
            n.object.position.addScaledVector(T, qe), n.object.updateMatrixWorld();
          } else if (n.object.isOrthographicCamera) {
            const Xe = new U(C.x, C.y, 0);
            Xe.unproject(n.object), n.object.zoom = Math.max(n.minZoom, Math.min(n.maxZoom, n.object.zoom / c)), n.object.updateProjectionMatrix(), de = !0;
            const qe = new U(C.x, C.y, 0);
            qe.unproject(n.object), n.object.position.sub(qe).add(Xe), n.object.updateMatrixWorld(), Fe = R.length();
          } else
            console.warn("WARNING: OrbitControls.js encountered an unknown camera type - zoom to cursor disabled."), n.zoomToCursor = !1;
          Fe !== null && (this.screenSpacePanning ? n.target.set(0, 0, -1).transformDirection(n.object.matrix).multiplyScalar(Fe).add(n.object.position) : (Qi.origin.copy(n.object.position), Qi.direction.set(0, 0, -1).transformDirection(n.object.matrix), Math.abs(n.object.up.dot(Qi.direction)) < Qp ? e.lookAt(n.target) : (Ka.setFromNormalAndCoplanarPoint(n.object.up, n.target), Qi.intersectPlane(Ka, n.target))));
        } else
          n.object.isOrthographicCamera && (n.object.zoom = Math.max(n.minZoom, Math.min(n.maxZoom, n.object.zoom / c)), n.object.updateProjectionMatrix(), de = !0);
        return c = 1, L = !1, de || he.distanceToSquared(n.object.position) > o || 8 * (1 - ge.dot(n.object.quaternion)) > o || De.distanceToSquared(n.target) > 0 ? (n.dispatchEvent(ja), he.copy(n.object.position), ge.copy(n.object.quaternion), De.copy(n.target), de = !1, !0) : !1;
      };
    }(), this.dispose = function() {
      n.domElement.removeEventListener("contextmenu", te), n.domElement.removeEventListener("pointerdown", xe), n.domElement.removeEventListener("pointercancel", Ce), n.domElement.removeEventListener("wheel", E), n.domElement.removeEventListener("pointermove", Me), n.domElement.removeEventListener("pointerup", Ce), n._domElementKeyEvents !== null && (n._domElementKeyEvents.removeEventListener("keydown", v), n._domElementKeyEvents = null);
    };
    const n = this, r = {
      NONE: -1,
      ROTATE: 0,
      DOLLY: 1,
      PAN: 2,
      TOUCH_ROTATE: 3,
      TOUCH_PAN: 4,
      TOUCH_DOLLY_PAN: 5,
      TOUCH_DOLLY_ROTATE: 6
    };
    let s = r.NONE;
    const o = 1e-6, a = new qa(), l = new qa();
    let c = 1;
    const h = new U(), f = new oe(), u = new oe(), m = new oe(), g = new oe(), x = new oe(), p = new oe(), d = new oe(), A = new oe(), _ = new oe(), T = new U(), C = new oe();
    let L = !1;
    const w = [], V = {};
    function M() {
      return 2 * Math.PI / 60 / 60 * n.autoRotateSpeed;
    }
    function y() {
      return Math.pow(0.95, n.zoomSpeed);
    }
    function Y(R) {
      l.theta -= R;
    }
    function ce(R) {
      l.phi -= R;
    }
    const B = function() {
      const R = new U();
      return function(_e, he) {
        R.setFromMatrixColumn(he, 0), R.multiplyScalar(-_e), h.add(R);
      };
    }(), H = function() {
      const R = new U();
      return function(_e, he) {
        n.screenSpacePanning === !0 ? R.setFromMatrixColumn(he, 1) : (R.setFromMatrixColumn(he, 0), R.crossVectors(n.object.up, R)), R.multiplyScalar(_e), h.add(R);
      };
    }(), G = function() {
      const R = new U();
      return function(_e, he) {
        const ge = n.domElement;
        if (n.object.isPerspectiveCamera) {
          const De = n.object.position;
          R.copy(De).sub(n.target);
          let Ve = R.length();
          Ve *= Math.tan(n.object.fov / 2 * Math.PI / 180), B(2 * _e * Ve / ge.clientHeight, n.object.matrix), H(2 * he * Ve / ge.clientHeight, n.object.matrix);
        } else
          n.object.isOrthographicCamera ? (B(_e * (n.object.right - n.object.left) / n.object.zoom / ge.clientWidth, n.object.matrix), H(he * (n.object.top - n.object.bottom) / n.object.zoom / ge.clientHeight, n.object.matrix)) : (console.warn("WARNING: OrbitControls.js encountered an unknown camera type - pan disabled."), n.enablePan = !1);
      };
    }();
    function Q(R) {
      n.object.isPerspectiveCamera || n.object.isOrthographicCamera ? c /= R : (console.warn("WARNING: OrbitControls.js encountered an unknown camera type - dolly/zoom disabled."), n.enableZoom = !1);
    }
    function X(R) {
      n.object.isPerspectiveCamera || n.object.isOrthographicCamera ? c *= R : (console.warn("WARNING: OrbitControls.js encountered an unknown camera type - dolly/zoom disabled."), n.enableZoom = !1);
    }
    function j(R) {
      if (!n.zoomToCursor)
        return;
      L = !0;
      const K = n.domElement.getBoundingClientRect(), _e = R.clientX - K.left, he = R.clientY - K.top, ge = K.width, De = K.height;
      C.x = _e / ge * 2 - 1, C.y = -(he / De) * 2 + 1, T.set(C.x, C.y, 1).unproject(e).sub(e.position).normalize();
    }
    function J(R) {
      return Math.max(n.minDistance, Math.min(n.maxDistance, R));
    }
    function ee(R) {
      f.set(R.clientX, R.clientY);
    }
    function I(R) {
      j(R), d.set(R.clientX, R.clientY);
    }
    function q(R) {
      g.set(R.clientX, R.clientY);
    }
    function pe(R) {
      u.set(R.clientX, R.clientY), m.subVectors(u, f).multiplyScalar(n.rotateSpeed);
      const K = n.domElement;
      Y(2 * Math.PI * m.x / K.clientHeight), ce(2 * Math.PI * m.y / K.clientHeight), f.copy(u), n.update();
    }
    function me(R) {
      A.set(R.clientX, R.clientY), _.subVectors(A, d), _.y > 0 ? Q(y()) : _.y < 0 && X(y()), d.copy(A), n.update();
    }
    function ve(R) {
      x.set(R.clientX, R.clientY), p.subVectors(x, g).multiplyScalar(n.panSpeed), G(p.x, p.y), g.copy(x), n.update();
    }
    function Ae(R) {
      j(R), R.deltaY < 0 ? X(y()) : R.deltaY > 0 && Q(y()), n.update();
    }
    function be(R) {
      let K = !1;
      switch (R.code) {
        case n.keys.UP:
          R.ctrlKey || R.metaKey || R.shiftKey ? ce(2 * Math.PI * n.rotateSpeed / n.domElement.clientHeight) : G(0, n.keyPanSpeed), K = !0;
          break;
        case n.keys.BOTTOM:
          R.ctrlKey || R.metaKey || R.shiftKey ? ce(-2 * Math.PI * n.rotateSpeed / n.domElement.clientHeight) : G(0, -n.keyPanSpeed), K = !0;
          break;
        case n.keys.LEFT:
          R.ctrlKey || R.metaKey || R.shiftKey ? Y(2 * Math.PI * n.rotateSpeed / n.domElement.clientHeight) : G(n.keyPanSpeed, 0), K = !0;
          break;
        case n.keys.RIGHT:
          R.ctrlKey || R.metaKey || R.shiftKey ? Y(-2 * Math.PI * n.rotateSpeed / n.domElement.clientHeight) : G(-n.keyPanSpeed, 0), K = !0;
          break;
      }
      K && (R.preventDefault(), n.update());
    }
    function Te() {
      if (w.length === 1)
        f.set(w[0].pageX, w[0].pageY);
      else {
        const R = 0.5 * (w[0].pageX + w[1].pageX), K = 0.5 * (w[0].pageY + w[1].pageY);
        f.set(R, K);
      }
    }
    function ke() {
      if (w.length === 1)
        g.set(w[0].pageX, w[0].pageY);
      else {
        const R = 0.5 * (w[0].pageX + w[1].pageX), K = 0.5 * (w[0].pageY + w[1].pageY);
        g.set(R, K);
      }
    }
    function Ye() {
      const R = w[0].pageX - w[1].pageX, K = w[0].pageY - w[1].pageY, _e = Math.sqrt(R * R + K * K);
      d.set(0, _e);
    }
    function we() {
      n.enableZoom && Ye(), n.enablePan && ke();
    }
    function b() {
      n.enableZoom && Ye(), n.enableRotate && Te();
    }
    function le(R) {
      if (w.length == 1)
        u.set(R.pageX, R.pageY);
      else {
        const _e = z(R), he = 0.5 * (R.pageX + _e.x), ge = 0.5 * (R.pageY + _e.y);
        u.set(he, ge);
      }
      m.subVectors(u, f).multiplyScalar(n.rotateSpeed);
      const K = n.domElement;
      Y(2 * Math.PI * m.x / K.clientHeight), ce(2 * Math.PI * m.y / K.clientHeight), f.copy(u);
    }
    function Z(R) {
      if (w.length === 1)
        x.set(R.pageX, R.pageY);
      else {
        const K = z(R), _e = 0.5 * (R.pageX + K.x), he = 0.5 * (R.pageY + K.y);
        x.set(_e, he);
      }
      p.subVectors(x, g).multiplyScalar(n.panSpeed), G(p.x, p.y), g.copy(x);
    }
    function re(R) {
      const K = z(R), _e = R.pageX - K.x, he = R.pageY - K.y, ge = Math.sqrt(_e * _e + he * he);
      A.set(0, ge), _.set(0, Math.pow(A.y / d.y, n.zoomSpeed)), Q(_.y), d.copy(A);
    }
    function $(R) {
      n.enableZoom && re(R), n.enablePan && Z(R);
    }
    function ye(R) {
      n.enableZoom && re(R), n.enableRotate && le(R);
    }
    function xe(R) {
      n.enabled !== !1 && (w.length === 0 && (n.domElement.setPointerCapture(R.pointerId), n.domElement.addEventListener("pointermove", Me), n.domElement.addEventListener("pointerup", Ce)), se(R), R.pointerType === "touch" ? O(R) : ze(R));
    }
    function Me(R) {
      n.enabled !== !1 && (R.pointerType === "touch" ? ie(R) : Ze(R));
    }
    function Ce(R) {
      Ee(R), w.length === 0 && (n.domElement.releasePointerCapture(R.pointerId), n.domElement.removeEventListener("pointermove", Me), n.domElement.removeEventListener("pointerup", Ce)), n.dispatchEvent(Za), s = r.NONE;
    }
    function ze(R) {
      let K;
      switch (R.button) {
        case 0:
          K = n.mouseButtons.LEFT;
          break;
        case 1:
          K = n.mouseButtons.MIDDLE;
          break;
        case 2:
          K = n.mouseButtons.RIGHT;
          break;
        default:
          K = -1;
      }
      switch (K) {
        case Pn.DOLLY:
          if (n.enableZoom === !1)
            return;
          I(R), s = r.DOLLY;
          break;
        case Pn.ROTATE:
          if (R.ctrlKey || R.metaKey || R.shiftKey) {
            if (n.enablePan === !1)
              return;
            q(R), s = r.PAN;
          } else {
            if (n.enableRotate === !1)
              return;
            ee(R), s = r.ROTATE;
          }
          break;
        case Pn.PAN:
          if (R.ctrlKey || R.metaKey || R.shiftKey) {
            if (n.enableRotate === !1)
              return;
            ee(R), s = r.ROTATE;
          } else {
            if (n.enablePan === !1)
              return;
            q(R), s = r.PAN;
          }
          break;
        default:
          s = r.NONE;
      }
      s !== r.NONE && n.dispatchEvent(Xr);
    }
    function Ze(R) {
      switch (s) {
        case r.ROTATE:
          if (n.enableRotate === !1)
            return;
          pe(R);
          break;
        case r.DOLLY:
          if (n.enableZoom === !1)
            return;
          me(R);
          break;
        case r.PAN:
          if (n.enablePan === !1)
            return;
          ve(R);
          break;
      }
    }
    function E(R) {
      n.enabled === !1 || n.enableZoom === !1 || s !== r.NONE || (R.preventDefault(), n.dispatchEvent(Xr), Ae(R), n.dispatchEvent(Za));
    }
    function v(R) {
      n.enabled === !1 || n.enablePan === !1 || be(R);
    }
    function O(R) {
      switch (ae(R), w.length) {
        case 1:
          switch (n.touches.ONE) {
            case Ln.ROTATE:
              if (n.enableRotate === !1)
                return;
              Te(), s = r.TOUCH_ROTATE;
              break;
            case Ln.PAN:
              if (n.enablePan === !1)
                return;
              ke(), s = r.TOUCH_PAN;
              break;
            default:
              s = r.NONE;
          }
          break;
        case 2:
          switch (n.touches.TWO) {
            case Ln.DOLLY_PAN:
              if (n.enableZoom === !1 && n.enablePan === !1)
                return;
              we(), s = r.TOUCH_DOLLY_PAN;
              break;
            case Ln.DOLLY_ROTATE:
              if (n.enableZoom === !1 && n.enableRotate === !1)
                return;
              b(), s = r.TOUCH_DOLLY_ROTATE;
              break;
            default:
              s = r.NONE;
          }
          break;
        default:
          s = r.NONE;
      }
      s !== r.NONE && n.dispatchEvent(Xr);
    }
    function ie(R) {
      switch (ae(R), s) {
        case r.TOUCH_ROTATE:
          if (n.enableRotate === !1)
            return;
          le(R), n.update();
          break;
        case r.TOUCH_PAN:
          if (n.enablePan === !1)
            return;
          Z(R), n.update();
          break;
        case r.TOUCH_DOLLY_PAN:
          if (n.enableZoom === !1 && n.enablePan === !1)
            return;
          $(R), n.update();
          break;
        case r.TOUCH_DOLLY_ROTATE:
          if (n.enableZoom === !1 && n.enableRotate === !1)
            return;
          ye(R), n.update();
          break;
        default:
          s = r.NONE;
      }
    }
    function te(R) {
      n.enabled !== !1 && R.preventDefault();
    }
    function se(R) {
      w.push(R);
    }
    function Ee(R) {
      delete V[R.pointerId];
      for (let K = 0; K < w.length; K++)
        if (w[K].pointerId == R.pointerId) {
          w.splice(K, 1);
          return;
        }
    }
    function ae(R) {
      let K = V[R.pointerId];
      K === void 0 && (K = new oe(), V[R.pointerId] = K), K.set(R.pageX, R.pageY);
    }
    function z(R) {
      const K = R.pointerId === w[0].pointerId ? w[1] : w[0];
      return V[K.pointerId];
    }
    n.domElement.addEventListener("contextmenu", te), n.domElement.addEventListener("pointerdown", xe), n.domElement.addEventListener("pointercancel", Ce), n.domElement.addEventListener("wheel", E, { passive: !1 }), this.update();
  }
}
const tm = "FeatureCollection", nm = [
  {
    type: "Feature",
    properties: {
      id_piso: 1207301801,
      cod_cat: 12073018,
      nivel: 1,
      altura: 2.5,
      area: 111.44
    },
    geometry: {
      type: "MultiPolygon",
      crs: {
        type: "name",
        properties: {
          name: "EPSG:4326"
        }
      },
      coordinates: [
        [
          [
            [
              -77.12030149593639,
              -11.987415412908778
            ],
            [
              -77.12035546688175,
              -11.987426229140224
            ],
            [
              -77.12035711928839,
              -11.987410951956447
            ],
            [
              -77.12037358088453,
              -11.987414250811023
            ],
            [
              -77.12037192847886,
              -11.987429527994742
            ],
            [
              -77.12038425169297,
              -11.987431997790118
            ],
            [
              -77.12045891273075,
              -11.987446960207757
            ],
            [
              -77.12046619084876,
              -11.98737965016659
            ],
            [
              -77.12030626370822,
              -11.9873647177975
            ],
            [
              -77.12030149593639,
              -11.987415412908778
            ]
          ]
        ]
      ]
    }
  },
  {
    type: "Feature",
    properties: {
      id_piso: 1207301802,
      cod_cat: 12073018,
      nivel: 2,
      altura: 2.5,
      area: 117.03
    },
    geometry: {
      type: "MultiPolygon",
      crs: {
        type: "name",
        properties: {
          name: "EPSG:4326"
        }
      },
      coordinates: [
        [
          [
            [
              -77.12030149607118,
              -11.987415413292887
            ],
            [
              -77.12035546700947,
              -11.987426228620599
            ],
            [
              -77.12035711941611,
              -11.987410951436821
            ],
            [
              -77.12037358009432,
              -11.987414250298345
            ],
            [
              -77.12037192861364,
              -11.98742952837885
            ],
            [
              -77.1203842518207,
              -11.987431997270495
            ],
            [
              -77.12045891286554,
              -11.987446960591848
            ],
            [
              -77.12046565877202,
              -11.987448312168949
            ],
            [
              -77.12047293688829,
              -11.98738100212789
            ],
            [
              -77.12046619097649,
              -11.987379649646982
            ],
            [
              -77.12030626384299,
              -11.987364718181604
            ],
            [
              -77.12030149607118,
              -11.987415413292887
            ]
          ]
        ]
      ]
    }
  }
], im = {
  type: tm,
  features: nm
};
var tn, _i, _n, rs, Go = 0;
const Ja = new $p(), yn = new oe(), Vo = ["#00a5e3", "#8dd7bf", "#6c88c4", "#ffa23a", "#ffd872", "#ff3d18", "#ff6446", "#ff8b74"];
class hm {
  constructor(e, t = "http://108.181.190.199/test/get_unidades.php?id=") {
    if (this.id = e, this.geojson = im, this.endpoint = t, e != 0) {
      this.geojson = rm(t, e).then((n) => {
        this.geojson = n, this.render();
      });
      return;
    }
    return this.render();
  }
  render() {
    _n = new Po({ antialias: !0 }), _n.shadowMap.enabled = !0, _n.setSize(window.innerWidth, window.innerHeight), document.body.appendChild(_n.domElement), tn = new mp(), _i = new At(45, window.innerWidth / window.innerHeight, 0.1, 1e3), rs = new em(_i, _n.domElement), this.x_values = [], this.z_values = [], this.axis1_min = 0, this.axis3_min = 0, this.axis1_max = 0, this.axis3_max = 0;
    var e = 5e3, t = 0.05;
    this.calculateMinAndMax(), _n.setClearColor(16777215, 0);
    const n = new Kp("white", 2);
    n.castShadow = !0, n.shadow.mapSize.width = 10240, n.shadow.mapSize.height = 10240, n.shadow.radius = 10;
    const r = new Jp(4210752, 20);
    tn.add(r), n.position.set(1, 2, 1), tn.add(n), this.geojson.features.forEach((l) => {
      l.geometry.coordinates.forEach((d) => {
        this.convertToPlaneCoords(d[0]);
      });
      let c = l.properties.nivel, h = l.properties.altura * t, u = this.getVertex(this.x_values, h, this.z_values).map((d) => new oe(d.x, d.z));
      const m = {
        steps: 1,
        depth: h,
        bevelEnabled: !1,
        bevelThickness: 1,
        bevelSize: 1,
        bevelOffset: 0,
        bevelSegments: 1
      }, g = new No(u), x = new gs(g, m), p = new Nt(
        x,
        new ns({ color: Vo[l.properties.nivel], side: Mt, transparent: !0, opacity: 0.8 })
      );
      p.scale.set(e, e, 1), p.rotateX(Math.PI / 2), p.position.set(0, h * c, 0), p.castShadow = !0, p.id_piso = l.properties.id_piso, p.nivel = c, tn.add(p), this.clearArrays();
    });
    const s = new ar(100, 100), o = new zo({ color: "white", side: Mt }), a = new Nt(s, o);
    a.receiveShadow = !0, a.rotateX(Math.PI / 2), tn.add(a), _i.position.set(0.5, 1, -1.5), rs.update(), ko();
  }
  convertToPlaneCoords(e) {
    e.forEach((t) => {
      let n = t[0], r = t[1];
      this.x_values.push(r), this.z_values.push(n);
    });
  }
  calculateMinAndMax() {
    let e = [], t = [], n = [];
    this.geojson.features.forEach((r) => {
      r.properties.nivel, n.push(r.properties.nivel), r.geometry.coordinates.forEach((s) => {
        s[0].forEach((o) => {
          e.push(o[1]), t.push(o[0]);
        });
      });
    }), this.axis1_min = Math.min(...e), this.axis3_min = Math.min(...t), this.axis1_max = Math.max(...e), this.axis3_max = Math.max(...t);
  }
  getVertex(e, t, n) {
    let r = [];
    for (var s = 0; s < e.length; s++) {
      let o = e[s] - (this.axis1_min + this.axis1_max) / 2, a = n[s] - (this.axis3_min + this.axis3_max) / 2;
      r.push(new U(o, t, a));
    }
    return r;
  }
  clearArrays() {
    this.x_values = [], this.z_values = [];
  }
}
async function rm(i, e) {
  try {
    return await (await fetch(i + e)).json();
  } catch {
    console.log("Error al obtener el GeoJSON desde el endpoint, verifique que el endpoint sea correcto");
  }
}
function ko() {
  requestAnimationFrame(ko), rs.update(), sm();
}
function sm() {
  om(), am(), _n.render(tn, _i);
}
function am() {
  tn.children.forEach((i) => {
    i.nivel != null && (i.id == Go ? i.material = new ns({ color: "red", side: Mt, transparent: !0, opacity: 0.8 }) : i.material = new ns({ color: Vo[i.nivel], side: Mt, transparent: !0, opacity: 0.8 }));
  });
}
function om() {
  if (yn.x != 0 && yn.y != 0) {
    Ja.setFromCamera(yn, _i);
    let i = Ja.intersectObjects(tn.children);
    i.length > 0 && i[0].object.id_piso != null && (console.log("Nivel: " + i[0].object.nivel + " ID: " + i[0].object.id_piso), Go = i[0].object.id);
  }
}
function lm(i) {
  yn.x = i.clientX / window.innerWidth * 2 - 1, yn.y = -(i.clientY / window.innerHeight) * 2 + 1;
}
function cm() {
  setTimeout(() => {
    yn.x = 0, yn.y = 0;
  }, 10);
}
window.addEventListener("mousedown", lm);
window.addEventListener("mouseup", cm);
export {
  hm as default
};
