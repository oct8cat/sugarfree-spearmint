export type Vec = [number, number, number];

export function cross(a: Vec, b: Vec): Vec {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

export function dot(a: Vec, b: Vec): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function mag([x, y, z]: Vec): number {
  return Math.sqrt(x ** 2 + y ** 2 + z ** 2);
}

export function scale([x, y, z]: Vec, s: number): Vec {
  return [x * s, y * s, z * s];
}

export function sub([ax, ay, az]: Vec, [bx, by, bz]: Vec): Vec {
  return [ax - bx, ay - by, az - bz];
}
