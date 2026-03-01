import { cross, dot, mag, scale, sub, type Vec } from "./vec";

export type Keplerian = {
  semiMajorAxis: number;
  eccentricity: number;
  inclination: number;
  raan: number;
  argumentOfPeriapsis: number;
  trueAnomaly: number;
};

export type State = {
  position: Vec;
  velocity: Vec;
};

export const TAU = Math.PI * 2;

/** Wraps an angle to [0, 2π] */
function wrap(angle: number) {
  return angle < 0 ? angle + TAU : angle % TAU;
}

/**
 * Converts state vectors to Keplerian elements.
 *
 * @param mu Standard gravitational parameter of the central body.
 * @param state Cartesian position and velocity.
 * @param k Z-axis unit vector
 * @param eps Error tolerance
 */
export function stateToKeplerian(mu: number, state: State, k: Vec = [0, 0, 1], eps: number = 1e-10): Keplerian {
  const { position: r, velocity: v } = state;
  const rmag = mag(r);
  const vmag = mag(v);
  const h = cross(r, v);
  const hmag = mag(h);
  const n = cross(k, h);
  const nmag = mag(n);
  const e = scale(sub(scale(r, vmag ** 2 - mu / rmag), scale(v, dot(r, v))), 1 / mu);
  const emag = mag(e);
  const a = 1 / (2 / rmag - vmag ** 2 / mu);
  const i = Math.atan2(Math.sqrt(h[0] ** 2 + h[1] ** 2), h[2]);
  const raan = nmag > eps ? wrap(Math.atan2(n[1], n[0])) : 0;
  const aop = nmag > eps && emag > eps ? wrap(Math.atan2(dot(cross(h, n), e), hmag * dot(n, e))) : 0;
  const ta = emag > eps ? wrap(Math.atan2((hmag * dot(r, v)) / (rmag * emag * mu), dot(e, r) / (emag * rmag))) : 0;

  return {
    semiMajorAxis: a,
    eccentricity: emag,
    inclination: i,
    raan,
    argumentOfPeriapsis: aop,
    trueAnomaly: ta,
  };
}

/**
 * Converts Keplerian elements to state vectors
 *
 * @param mu Standard gravitational parameter of the central body.
 * @param keplerian Keplerian elements to convert
 */

export function keplerianToState(mu: number, keplerian: Keplerian): State {
  const {
    semiMajorAxis: a,
    eccentricity: e,
    inclination: i,
    raan,
    argumentOfPeriapsis: aop,
    trueAnomaly: ta,
  } = keplerian;

  const cosTa = Math.cos(ta);
  const sinTa = Math.sin(ta);
  const cosRaan = Math.cos(raan);
  const sinRaan = Math.sin(raan);
  const cosI = Math.cos(i);
  const sinI = Math.sin(i);
  const cosAop = Math.cos(aop);
  const sinAop = Math.sin(aop);
  const p = Math.abs(a * (1 - e ** 2));
  const rmag = p / (1 + e * cosTa);
  const x = rmag * cosTa;
  const y = rmag * sinTa;
  const vscale = Math.sqrt(mu / Math.max(p, 1e-12));
  const vx = vscale * -sinTa;
  const vy = vscale * (e + cosTa);

  const r: Vec = [
    x * (cosRaan * cosAop - sinRaan * cosI * sinAop) - y * (cosRaan * sinAop + sinRaan * cosI * cosAop),
    x * (sinRaan * cosAop + cosRaan * cosI * sinAop) - y * (sinRaan * sinAop - cosRaan * cosI * cosAop),
    x * (sinI * sinAop) + y * (sinI * cosAop),
  ];
  const v: Vec = [
    vx * (cosRaan * cosAop - sinRaan * cosI * sinAop) - vy * (cosRaan * sinAop + sinRaan * cosI * cosAop),
    vx * (sinRaan * cosAop + cosRaan * cosI * sinAop) - vy * (sinRaan * sinAop - cosRaan * cosI * cosAop),
    vx * (sinI * sinAop) + vy * (sinI * cosAop),
  ];

  return { position: r, velocity: v };
}
