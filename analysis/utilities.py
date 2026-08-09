from nt import truncate
from numpy import *
import numpy as np

from scipy.signal import savgol_filter

def unwrap_phase(phi):
    """Supposed monotonic increasing phase.
    """
    diff_phase = diff(phi)
    unwPhase = cumsum(r_[0, abs(diff_phase)]) #Add the 0 

    return unwPhase

# def unwrap_phase(phi):
#     phi = (phi ) % (pi) - pi/2
#     return phi
#     d = diff(phi)
#     #d = (d ) % (pi)
#     return cumsum(r_[phi[0], d])

def unwrap_increasing(phi):
    d = diff(phi)
    d = (d + pi) % (2*pi) - pi
    d = maximum(d, 0.0)          # ignore small backward noise
    return cumsum(r_[phi[0], d])
    
def norm_I(I, filter = True, polyorder = 5):
    if filter:
        I = savgol_filter(I, window_length=11, polyorder=polyorder)

    I_norm = 2 * (I - I.min()) / (I.max() - I.min()) - 1
    return I_norm

def phase_unwrap_from_I(I, filter = True, polyorder = 5):
    
    if filter:
        I_norm = savgol_filter(I, window_length=11, polyorder=polyorder)
        
    I_norm = norm_I(I)
    phase = arccos(I_norm)
    return unwrap_phase(phase)


def fit_ellipse_conic(x, y):
    """Direct least-squares ellipse (Fitzgibbon / Halir).

    Fits  A x^2 + B xy + C y^2 + D x + E y + F = 0
    with constraint  4AC - B^2 = 1.

    Returns
    -------
    A, B, C, D, E, F : float
    """
    x = asarray(x, dtype=float).ravel()
    y = asarray(y, dtype=float).ravel()
    if x.size != y.size:
        raise ValueError('x and y must have the same length')
    if x.size < 6:
        raise ValueError('need at least 6 points')

    D1 = c_[x**2, x*y, y**2]
    D2 = c_[x, y, ones_like(x)]
    S1 = D1.T @ D1
    S2 = D1.T @ D2
    S3 = D2.T @ D2
    T = -linalg.inv(S3) @ S2.T
    M = S1 + S2 @ T
    M = array([M[2, :] / 2.0,
               -M[1, :],
               M[0, :] / 2.0])

    _, evecs = linalg.eig(M)
    cond = 4.0 * evecs[0, :] * evecs[2, :] - evecs[1, :]**2
    i_ell = where(cond > 0)[0]
    if i_ell.size == 0:
        raise RuntimeError('No ellipse solution (hyperbola/parabola). Plot Ix vs Iy.')

    a1 = evecs[:, i_ell[0]].real
    A, B, C, D, E, F = r_[a1, (T @ a1).ravel()].real
    return A, B, C, D, E, F

def conic_to_heydemann(A, B, C, D, E, F):
    """Conic ABCDEF -> Heydemann (x0, y0, a, b, eps).

    Ix = x0 + a cos(phi)
    Iy = y0 + b cos(phi - delta),   eps = cos(delta)

    a = amplitude of Ix,  b = amplitude of Iy.
    """
    # positive-definite quadratic form (A>0, C>0)
    if A < 0:
        A, B, C, D, E, F = -A, -B, -C, -D, -E, -F

    den = B**2 - 4.0*A*C          # < 0 for ellipse
    if den >= 0:
        raise ValueError('not an ellipse (B^2-4AC >= 0)')

    x0 = (2.0*C*D - B*E) / den
    y0 = (2.0*A*E - B*D) / den

    Fp = A*x0**2 + B*x0*y0 + C*y0**2 + D*x0 + E*y0 + F
    if Fp >= 0:
        raise ValueError('F\' >= 0 after sign flip; check the conic fit')

    eps = -B / (2.0 * sqrt(A*C))
    eps = clip(eps, -1.0 + 1e-12, 1.0 - 1e-12)
    one_e2 = 1.0 - eps**2

    a = sqrt((-Fp) / (A * one_e2))   # Ix
    b = sqrt((-Fp) / (C * one_e2))   # Iy
    return x0, y0, a, b, eps

def conic_to_geometry(A, B, C, D, E, F):
    """Geometric ellipse: centre, major-axis angle (rad), semi-major, semi-minor."""
    if A < 0:
        A, B, C, D, E, F = -A, -B, -C, -D, -E, -F
    den = B**2 - 4.0*A*C
    x0 = (2.0*C*D - B*E) / den
    y0 = (2.0*A*E - B*D) / den
    Fp = A*x0**2 + B*x0*y0 + C*y0**2 + D*x0 + E*y0 + F

    theta = 0.5 * arctan2(B, A - C)          # x'-axis after removing xy term
    ct, st = cos(theta), sin(theta)
    Ap = A*ct**2 + C*st**2 + B*st*ct
    Cp = A*st**2 + C*ct**2 - B*st*ct
    rx = sqrt((-Fp) / Ap)                     # along theta
    ry = sqrt((-Fp) / Cp)                     # along theta + pi/2

    if rx >= ry:
        major, minor, theta_maj = rx, ry, theta
    else:
        major, minor, theta_maj = ry, rx, theta + 0.5*pi
    theta_maj = (theta_maj + 0.5*pi) % pi - 0.5*pi   # (-90, 90]
    return x0, y0, theta_maj, major, minor


def phase_from_heydemann(Ix, Iy, x0, y0, a, b, eps):
    """Invert Heydemann ellipse to wrapped phi in (-pi, pi].

    Matches conic_to_heydemann:
        Ix = x0 + a cos(phi)
        Iy = y0 + b [eps cos(phi) + sqrt(1-eps^2) sin(phi)]

    sin_sign: +1 or -1 if phi(V) comes out decreasing.
    """
    X = (asarray(Ix, float).ravel() - x0) / a
    Y = (asarray(Iy, float).ravel() - y0) / b
    s = sqrt(1.0 - eps**2)
    cphi = X
    sphi = 1 * (Y - eps * X) / s

    phi_v = arctan2(sphi, cphi)
    phi_v_unwrap = unwrap(phi_v)
    #If the phase is decreasing with V, we need to invert the sign of the sin term - Recalculate
    if phi_v_unwrap[-1] < phi_v_unwrap[0]:         
        sphi = -1 * (Y - eps * X) / s    # SLM phase should increase with V
        phi_v = arctan2(sphi, cphi)

    return phi_v

def visualize_conic_fit(Ix, Iy, A, B, C, D, E, F, x0=None, y0=None,
                        a=None, b=None, eps=None, ax=None):
    """Lissajous + fitted ellipse, centre, major/minor, theta, epsilon."""
    import matplotlib.pyplot as plt

    Ix = asarray(Ix, float).ravel()
    Iy = asarray(Iy, float).ravel()

    x0_g, y0_g, theta_maj, major, minor = conic_to_geometry(A, B, C, D, E, F)
    if x0 is None:
        x0 = x0_g
    if y0 is None:
        y0 = y0_g
    if a is None or b is None or eps is None:
        _, _, a, b, eps = conic_to_heydemann(A, B, C, D, E, F)

    pad = 0.1 * maximum(ptp(Ix), ptp(Iy))
    xx = linspace(Ix.min() - pad, Ix.max() + pad, 400)
    yy = linspace(Iy.min() - pad, Iy.max() + pad, 400)
    XX, YY = meshgrid(xx, yy)
    Q = A*XX**2 + B*XX*YY + C*YY**2 + D*XX + E*YY + F

    ct, st = cos(theta_maj), sin(theta_maj)
    if ax is None:
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
    else:
        fig = ax.figure

    ax.plot(Ix, Iy, '.', ms=4, label='data')
    ax.contour(XX, YY, Q, levels=[0], colors='r', linewidths=2)
    ax.plot([x0 - major*ct, x0 + major*ct],
            [y0 - major*st, y0 + major*st], 'g-', lw=2, label='major')
    ax.plot([x0 + minor*st, x0 - minor*st],
            [y0 - minor*ct, y0 + minor*ct], 'b-', lw=2, label='minor')
    ax.plot(x0, y0, 'ko', ms=8, label='centre')

    r_ang = 0.35 * major
    ax.plot([x0, x0 + r_ang], [y0, y0], 'k--', lw=1)
    ax.plot([x0, x0 + r_ang*ct], [y0, y0 + r_ang*st], 'k--', lw=1)
    ax.text(x0 + 0.15*major, y0 + 0.08*major,
            rf'$\theta={theta_maj*180/pi:.1f}^\circ$', fontsize=11)

    ax.annotate(
        rf'centre$=({x0:.3f},\ {y0:.3f})$',
        xy=(x0, y0), xytext=(12, -18), textcoords='offset points',
        fontsize=10, arrowprops=dict(arrowstyle='->', color='k', lw=0.8),
    )
    ax.text(0.03, 0.97,
            rf'$\varepsilon={eps:.3f}$' + '\n' + rf'$a={a:.3f},\ b={b:.3f}$',
            transform=ax.transAxes, va='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='w', alpha=0.8))

    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel(r'$I_x$ (H)')
    ax.set_ylabel(r'$I_y$ (V)')
    ax.set_title('Conic / ellipse fit')
    ax.legend(loc='best')
    fig.tight_layout()
    return ax


def quadrature_detection(Ix, Iy, ax = None, debugging = True):
    A, B, C, D, E, F = fit_ellipse_conic(Ix, Iy)
    x0, y0, a, b, eps = conic_to_heydemann(A, B, C, D, E, F)

    if debugging:
        print(f'x0={x0:.4g}  y0={y0:.4g}  a={a:.4g}  b={b:.4g}  eps={eps:.4g}')
        disc = B**2 - 4*A*C
        print(f'A={A:.6g}  B={B:.6g}  C={C:.6g}  D={D:.6g}  E={E:.6g}  F={F:.6g}')
        print(f'B^2 - 4AC = {disc:.6g}   (must be < 0 for ellipse)')

        visualize_conic_fit(Ix, Iy, A, B, C, D, E, F, x0=x0, y0=y0, ax = ax)
    
    phi = phase_from_heydemann(Ix, Iy, x0, y0, a, b, eps)
    phi_unw = unwrap(phi)   
    phi_v_Quad =  phi_unw - phi_unw[0]

    return phi_v_Quad


    