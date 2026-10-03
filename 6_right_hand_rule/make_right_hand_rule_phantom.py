"""Warp the healthy-adult hand template into a right-hand-rule pose.

Source: Hegdé et al., the NIfTI behind NIH 3D entry 17237
(https://3d.nih.gov/entries/17237), distributed from
https://github.com/HegdeUSA/Hand_template (MIT).

Target pose, matching the usual drawing of the rule (palm in the plane of
the page, index up, thumb to the right, the toward-the-viewer axis on the
middle finger):

- palm, index finger, and thumb lie in the world XY plane
- index finger points along +Y
- thumb points along +X
- middle finger points along +Z

World +X, +Y, +Z are a right-handed triad: thumb × index = middle finger.
Ring and pinky are flexed toward the palm, as in that drawing.
"""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

# Labels in the source grid. Track ids during the walk are separate from these.
PALM = 1
THUMB = 2
INDEX = 3
MIDDLE = 4
RING = 5
PINKY = 6
_BROAD = np.uint8(255)
FINGER_NAMES = {
    THUMB: "thumb",
    INDEX: "index",
    MIDDLE: "middle",
    RING: "ring",
    PINKY: "pinky",
}


def load_hand(path, step=2, threshold=500.0):
    """Crop the template and subsample to about 1 mm.

    The file is 0.5 mm and already aligned with the array axes: k runs
    wrist → fingertips, i runs pinky → thumb, j is the palmar/dorsal thickness.
    """
    img = nib.load(str(path))
    data = np.asanyarray(img.dataobj, dtype=np.float32)
    zooms = np.asarray(img.header.get_zooms()[:3], dtype=float)
    mask_full = data > threshold
    coords = np.where(mask_full)
    bb = tuple(slice(int(a.min()), int(a.max()) + 1) for a in coords)
    data = data[bb][::step, ::step, ::step]
    spacing = zooms * float(step)
    # The template axes are within a degree of orthogonal. Use the mean
    # spacing so the phantom is isotropic in its own world frame.
    voxel_mm = float(np.mean(spacing))
    mask = data > threshold
    mask = ndimage.binary_fill_holes(mask)
    mask = ndimage.binary_opening(mask, iterations=1)
    labeled, n = ndimage.label(mask)
    if n > 1:
        sizes = ndimage.sum(mask, labeled, index=range(1, n + 1))
        mask = labeled == (int(np.argmax(sizes)) + 1)
    return data, mask, voxel_mm


def _components(slice_2d, min_size):
    lab, n = ndimage.label(slice_2d)
    comps = []
    for c in range(1, n + 1):
        ii, jj = np.where(lab == c)
        if ii.size < min_size:
            continue
        comps.append((float(ii.mean()), ii, jj))
    comps.sort(key=lambda item: item[0])
    return comps


def label_digits(mask):
    """Separate thumb and the four fingers from the palm.

    Walking from the fingertips back to the wrist, a small cross-section is a
    digit and a large one is palm. The thumb is the radial digit that never
    reaches the fingertips.
    """
    i_size, j_size, k_size = mask.shape
    labels = np.zeros(mask.shape, np.uint8)
    # Finger cross-sections at 1 mm are a few hundred voxels. The palm is thousands.
    finger_max = 700
    match_dist = 16.0
    min_size = 25

    active = []  # tracks still being extended
    finished = []
    next_id = 1

    def new_track(mean_i, k):
        nonlocal next_id
        track = {
            "id": next_id,
            "last_i": mean_i,
            "kmax": k,
            "kmin": k,
            "sum_i": 0.0,
            "count": 0,
        }
        next_id += 1
        return track

    for k in range(k_size - 1, -1, -1):
        comps = _components(mask[:, :, k], min_size)
        used = set()
        still = []
        for mean_i, ii, jj in comps:
            if ii.size > finger_max:
                labels[ii, jj, k] = _BROAD
                continue
            best = None
            best_d = match_dist
            for ai, track in enumerate(active):
                if ai in used:
                    continue
                dist = abs(track["last_i"] - mean_i)
                if dist < best_d:
                    best_d = dist
                    best = ai
            if best is None:
                track = new_track(mean_i, k)
            else:
                used.add(best)
                track = active[best]
                track["last_i"] = mean_i
                track["kmin"] = k
            labels[ii, jj, k] = track["id"]
            track["sum_i"] += mean_i * ii.size
            track["count"] += ii.size
            still.append(track)
        for ai, track in enumerate(active):
            if ai not in used:
                finished.append(track)
        active = still
    finished.extend(active)

    tracks = [t for t in finished if t["count"] > 400 and (t["kmax"] - t["kmin"]) > 8]
    if not tracks:
        raise RuntimeError("no digit tracks found in the hand template")

    # The four fingers reach well past the thumb. The pinky is the short one,
    # so this cut is the knuckle line rather than "near the longest tip".
    fingers = [t for t in tracks if t["kmax"] >= 160 and t["count"] > 600]
    fingers.sort(key=lambda t: t["sum_i"] / t["count"])
    if len(fingers) != 4:
        summary = [
            (round(t["sum_i"] / t["count"], 1), t["kmin"], t["kmax"], t["count"])
            for t in tracks
        ]
        raise RuntimeError(f"expected 4 fingers, found {len(fingers)} from {summary}")

    # Low i is ulnar (pinky). High i is radial (index, then the thumb).
    order = [PINKY, RING, MIDDLE, INDEX]
    remap = {track["id"]: name for track, name in zip(fingers, order)}
    finger_ids = set(remap)
    mean_i = {t["id"]: t["sum_i"] / t["count"] for t in tracks}
    index_i = mean_i[fingers[-1]["id"]]
    thumb_tracks = [
        t
        for t in tracks
        if t["id"] not in finger_ids
        and mean_i[t["id"]] > index_i + 8
        and t["kmax"] < 160
    ]
    if not thumb_tracks:
        raise RuntimeError("thumb track not found")
    thumb = max(thumb_tracks, key=lambda t: t["count"])
    remap[thumb["id"]] = THUMB

    out = np.zeros(mask.shape, np.uint8)
    for old, new in remap.items():
        out[labels == old] = new
    named = np.isin(labels, list(remap.keys()))
    out[(labels == _BROAD) | ((labels == 0) & mask) | ((labels > 0) & ~named)] = PALM
    return out


def _palmar_is_high_j(labels):
    """True when source +j points out of the palm.

    The four fingers lie on one face of the hand. The thick pad (thenar and
    the rest of the palm) is the other face, and that face is the front.
    Bright-voxel tests are not used: on this T1 template fat is as bright as
    marrow, and that test put the dorsum on the palmar side.
    """
    fingers = np.isin(labels, (INDEX, MIDDLE, RING, PINKY))
    finger_j = float(np.where(fingers)[1].mean())
    palm_j = float(np.where(labels == PALM)[1].mean())
    return palm_j > finger_j


def _normalize(v):
    n = np.linalg.norm(v)
    if n < 1e-8:
        raise RuntimeError("zero-length anatomical axis")
    return v / n


def _coords(labels, code, voxel_mm):
    ijk = np.stack(np.where(labels == code), axis=1).astype(np.float64)
    return ijk * voxel_mm


def _tip_and_base(points, distal):
    along = points @ distal
    base = points[along <= np.percentile(along, 8)].mean(axis=0)
    tip = points[along >= np.percentile(along, 98)].mean(axis=0)
    return tip, base


def _rotate_about_z(rel, angles):
    c = np.cos(angles)
    s = np.sin(angles)
    x = rel[:, 0] * c - rel[:, 1] * s
    y = rel[:, 0] * s + rel[:, 1] * c
    out = np.empty_like(rel)
    out[:, 0] = x
    out[:, 1] = y
    out[:, 2] = rel[:, 2]
    return out


def _blend(distance, length):
    t = np.clip(distance / length, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _remove_bow(rel, along, axis):
    """Drop a digit's longitudinal bend along `along` for one axis.

    Thickness stays. The template is not perfectly flat; taking the curve out
    keeps a straight digit in the plane it is supposed to occupy.
    """
    span = max(float(along.max() - along.min()), 1.0)
    nb = 48
    bins = np.clip(((along - along.min()) / span * (nb - 1)).astype(np.int32), 0, nb - 1)
    summed = np.zeros(nb, np.float64)
    cnt = np.zeros(nb, np.float64)
    np.add.at(summed, bins, rel[:, axis])
    np.add.at(cnt, bins, 1.0)
    known = cnt > 0
    idx = np.arange(nb)
    mean = np.interp(idx, idx[known], summed[known] / cnt[known])
    mean = ndimage.gaussian_filter1d(mean, 1.5)
    out = rel.copy()
    out[:, axis] = rel[:, axis] - mean[bins]
    return out


def _sweep_toward_palm(rel, theta_at):
    """Lay a digit onto a smooth arc in the YZ plane.

    `theta_at(y)` is the flexion angle in radians. Zero keeps the digit along
    +Y; a positive angle turns it toward +Z, the front of the hand. Each
    cross-section follows the arc, so the finger stays one piece.
    """
    flat = _remove_bow(rel, rel[:, 1], 2)
    y = flat[:, 1]
    y_min = float(y.min())
    y_max = float(y.max())
    grid = np.linspace(y_min, y_max, 800)
    # Proximal to the knuckle the tangent stays +Y, so the digit meets the palm.
    theta = theta_at(np.maximum(grid, 0.0))
    step = np.diff(grid, prepend=grid[0])
    cy = np.cumsum(np.cos(theta) * step)
    cz = np.cumsum(np.sin(theta) * step)
    cy -= np.interp(0.0, grid, cy)
    cz -= np.interp(0.0, grid, cz)
    cy_p = np.interp(y, grid, cy)
    cz_p = np.interp(y, grid, cz)
    theta_p = np.interp(y, grid, theta)
    oz = flat[:, 2]
    out = np.empty_like(flat)
    out[:, 0] = flat[:, 0]
    out[:, 1] = cy_p - oz * np.sin(theta_p)
    out[:, 2] = cz_p + oz * np.cos(theta_p)
    return out


def _aim_toward_center(rel, degrees, blend_length):
    """Yaw a straight digit about Z. Negative degrees carry the tip toward +X.

    +X is radial, toward the middle of the palm from the little finger.
    The angle grows away from the knuckle, so the base stays put.
    """
    along = rel[:, 1]
    ang = np.deg2rad(degrees) * _blend(along - along.min(), blend_length)
    c = np.cos(ang)
    s = np.sin(ang)
    out = rel.copy()
    out[:, 0] = rel[:, 0] * c - rel[:, 1] * s
    out[:, 1] = rel[:, 0] * s + rel[:, 1] * c
    return out


def _curl_angle(distance, length, mcp_deg, pip_deg):
    """Knuckle flexion, then a second bend so the tip comes back toward the palm."""
    mcp = np.deg2rad(mcp_deg) * _blend(distance, 0.38 * length)
    pip = np.deg2rad(pip_deg) * _blend(distance - 0.42 * length, 0.30 * length)
    return mcp + pip


def _pose_digits(local, codes):
    """Hinge each digit. `local` is Nx3 in the palm frame (mm).

    +X radial (thumb side), +Y distal, +Z palmar (front of the hand).
    Flexion is toward +Z. The index and the thumb stay in the XY plane.
    """
    posed = local.copy()

    def part(code):
        sel = codes == code
        return sel, posed[sel]

    # Index: swing it in the plane until it points along +Y, and take out the
    # dorsal bow so the finger lies on the front of the hand.
    sel, pts = part(INDEX)
    base = pts[pts[:, 1] <= np.percentile(pts[:, 1], 8)].mean(axis=0)
    tip = pts[pts[:, 1] >= np.percentile(pts[:, 1], 98)].mean(axis=0)
    direction = tip - base
    direction[2] = 0.0
    direction = _normalize(direction)
    # Tilt from +Y toward +X. A positive rotation about Z swings +X toward +Y.
    phi = np.arctan2(direction[0], direction[1])
    rel = pts - base
    rel = _remove_bow(rel, rel @ direction, 2)
    distal = rel[:, 1]
    angles = phi * _blend(distal - distal.min(), 18.0)
    posed[sel] = _rotate_about_z(rel, angles) + base

    # Thumb: same treatment, then swing it in the plane until it points along +X.
    sel, pts = part(THUMB)
    base = pts[pts[:, 0] <= np.percentile(pts[:, 0], 12)].mean(axis=0)
    tip = pts[np.linalg.norm(pts - base, axis=1).argmax()]
    direction = tip - base
    direction[2] = 0.0
    direction = _normalize(direction)
    phi = np.arctan2(direction[1], direction[0])  # angle from +X toward +Y
    rel = pts - base
    rel = _remove_bow(rel, rel @ direction, 2)
    along = rel @ direction
    angles = -phi * _blend(along - along.min(), 16.0)
    posed[sel] = _rotate_about_z(rel, angles) + base

    # Middle finger: 90° palmar flexion. Distal to the knuckle it is straight
    # along +Z, toward the viewer when the palm faces them.
    sel, pts = part(MIDDLE)
    base = pts[pts[:, 1] <= np.percentile(pts[:, 1], 8)].mean(axis=0)
    rel = pts - base
    posed[sel] = _sweep_toward_palm(rel, lambda d: (np.pi / 2.0) * _blend(d, 30.0)) + base

    # Ring folds down onto the palm, past the curl already on main.
    _curl_one(part, posed, RING, 112.0, 98.0, aim_deg=0.0)
    # Little finger folds with it and yaws further toward the middle of the palm.
    # More yaw than this runs the shaft through the ring finger.
    _curl_one(part, posed, PINKY, 112.0, 96.0, aim_deg=-28.0)

    return posed


def _curl_one(part, posed, code, mcp_deg, pip_deg, aim_deg):
    sel, pts = part(code)
    if pts.size == 0:
        return
    base = pts[pts[:, 1] <= np.percentile(pts[:, 1], 8)].mean(axis=0)
    rel = pts - base
    if aim_deg:
        rel = _aim_toward_center(rel, aim_deg, 34.0)
    length = max(float(np.percentile(rel[:, 1], 98)), 1.0)
    posed[sel] = _sweep_toward_palm(
        rel,
        lambda d, length=length, mcp_deg=mcp_deg, pip_deg=pip_deg: _curl_angle(
            d, length, mcp_deg, pip_deg
        ),
    ) + base


def _splat(points, values, origin, shape):
    acc = np.zeros(shape, np.float64)
    weight = np.zeros(shape, np.float64)
    q = points - origin
    i0 = np.floor(q).astype(np.int32)
    frac = q - i0
    for di in (0, 1):
        wx = frac[:, 0] if di else 1.0 - frac[:, 0]
        for dj in (0, 1):
            wy = frac[:, 1] if dj else 1.0 - frac[:, 1]
            for dk in (0, 1):
                wz = frac[:, 2] if dk else 1.0 - frac[:, 2]
                w = wx * wy * wz
                ii = i0[:, 0] + di
                jj = i0[:, 1] + dj
                kk = i0[:, 2] + dk
                valid = (
                    (w > 0)
                    & (ii >= 0)
                    & (jj >= 0)
                    & (kk >= 0)
                    & (ii < shape[0])
                    & (jj < shape[1])
                    & (kk < shape[2])
                )
                np.add.at(acc, (ii[valid], jj[valid], kk[valid]), values[valid] * w[valid])
                np.add.at(weight, (ii[valid], jj[valid], kk[valid]), w[valid])
    out = np.zeros(shape, np.float32)
    occupied = weight > 1e-6
    out[occupied] = (acc[occupied] / weight[occupied]).astype(np.float32)
    return out


def _drop_specks(volume, min_voxels=200):
    support = volume > 0
    labels, n = ndimage.label(support)
    if n <= 1:
        return volume
    sizes = ndimage.sum(support, labels, index=range(1, n + 1))
    keep = np.zeros(n + 1, dtype=bool)
    keep[1:] = np.asarray(sizes) >= min_voxels
    keep[int(np.argmax(sizes)) + 1] = True
    cleaned = volume.copy()
    cleaned[~keep[labels]] = 0
    return cleaned


def _fill_knuckle_gaps(volume):
    support = volume > 0
    closed = ndimage.binary_closing(support, structure=np.ones((3, 3, 3)))
    holes = closed & ~support
    if not np.any(holes):
        return volume
    _, inds = ndimage.distance_transform_edt(~support, return_indices=True)
    filled = volume.copy()
    filled[holes] = volume[inds[0][holes], inds[1][holes], inds[2][holes]]
    return filled


def _axis_report(posed, codes, distal_score):
    """Tip-minus-base vector of each posed digit, in mm.

    `distal_score` is the along-finger coordinate from before the hinge, so a
    curled fingertip is still the tip. The thumb is scored along +X instead.
    """
    report = {}
    for code, name in FINGER_NAMES.items():
        sel = codes == code
        pts = posed[sel]
        if pts.size == 0:
            continue
        score = pts[:, 0] if name == "thumb" else distal_score[sel]
        base = pts[score <= np.percentile(score, 8)].mean(axis=0)
        tip = pts[score >= np.percentile(score, 98)].mean(axis=0)
        report[name] = {
            "base_mm": base,
            "tip_mm": tip,
            "direction_mm": tip - base,
        }
    return report


def _check_pose(report):
    index = report["index"]["direction_mm"]
    thumb = report["thumb"]["direction_mm"]
    middle = report["middle"]["direction_mm"]
    if not (index[1] > 40 and abs(index[0]) < 12 and abs(index[2]) < 20):
        raise RuntimeError(f"index finger is not along +Y: {index}")
    if not (thumb[0] > 25 and abs(thumb[1]) < 30 and abs(thumb[2]) < 25):
        raise RuntimeError(f"thumb is not along +X: {thumb}")
    if not (middle[2] > 40 and abs(middle[0]) < 25 and abs(middle[1]) < 30):
        raise RuntimeError(f"middle finger is not along +Z: {middle}")
    # Ring and pinky curl back over the palm: the tip is proximal of the
    # knuckle and still on the palmar side of the hand.
    for name in ("ring", "pinky"):
        curled = report[name]["direction_mm"]
        if not (curled[2] > 5 and curled[1] < -14):
            raise RuntimeError(f"{name} is not curled into the palm: {curled}")
    # The little finger's tip moves toward +X, the middle of the palm.
    pinky = report["pinky"]["direction_mm"]
    if pinky[0] < 4:
        raise RuntimeError(f"pinky does not aim toward the center of the palm: {pinky}")
    # thumb × index should point the same way as the middle finger.
    cross = np.cross(thumb / np.linalg.norm(thumb), index / np.linalg.norm(index))
    if np.dot(cross, middle) <= 0:
        raise RuntimeError("posed digits are not a right-handed triad")


def build_phantom(nii_path, step=2, threshold=500.0):
    """Return (volume, affine, report) for the posed hand.

    `affine` maps voxel indices to millimetres with +X thumb, +Y index,
    +Z middle finger. The determinant is positive.
    """
    data, mask, voxel_mm = load_hand(nii_path, step=step, threshold=threshold)
    labels = label_digits(mask)
    labels[~mask] = 0

    palmar_is_high_j = _palmar_is_high_j(labels)
    palmar = np.array([0.0, 1.0 if palmar_is_high_j else -1.0, 0.0])

    middle_pts = _coords(labels, MIDDLE, voxel_mm)
    index_pts = _coords(labels, INDEX, voxel_mm)
    pinky_pts = _coords(labels, PINKY, voxel_mm)
    thumb_pts = _coords(labels, THUMB, voxel_mm)
    distal = _normalize(_tip_and_base(middle_pts, np.array([0.0, 0.0, 1.0]))[0]
                        - _tip_and_base(middle_pts, np.array([0.0, 0.0, 1.0]))[1])
    # Radio-ulnar axis, toward the thumb, lying in the palm.
    across = index_pts.mean(axis=0) - pinky_pts.mean(axis=0)
    across = across - distal * np.dot(across, distal)
    radial = _normalize(across)
    # Right-handed palmar normal from the template's own chirality.
    z_from_hand = _normalize(np.cross(radial, distal))
    # Keep the thumb on +X. If that right-handed normal points dorsally, mirror
    # Z so the palm faces +Z. The drawing is this view: palm toward the viewer,
    # thumb to the viewer's right, index up.
    mirror_z = bool(np.dot(z_from_hand, palmar) < 0)

    origin = _tip_and_base(middle_pts, distal)[1]
    basis = np.stack([radial, distal, z_from_hand], axis=1)  # columns are world axes

    ijk = np.stack(np.where(mask), axis=1).astype(np.float64)
    values = data[mask].astype(np.float32)
    codes = labels[mask]
    source_mm = ijk * voxel_mm
    local = (source_mm - origin) @ basis
    if mirror_z:
        local[:, 2] *= -1.0
    posed = _pose_digits(local, codes)
    report = _axis_report(posed, codes, local[:, 1])
    _check_pose(report)
    report["voxel_mm"] = voxel_mm
    report["mirrored_to_match_drawing"] = mirror_z
    report["palmar_is_high_j"] = bool(palmar_is_high_j)
    report["n_voxels"] = {
        "palm": int(np.sum(codes == PALM)),
        "thumb": int(np.sum(codes == THUMB)),
        "index": int(np.sum(codes == INDEX)),
        "middle": int(np.sum(codes == MIDDLE)),
        "ring": int(np.sum(codes == RING)),
        "pinky": int(np.sum(codes == PINKY)),
    }

    margin = 6.0
    lower = posed.min(axis=0) - margin
    upper = posed.max(axis=0) + margin
    shape = tuple(int(np.ceil(n)) for n in (upper - lower))
    volume = _splat(posed, values, lower, shape)
    volume = _fill_knuckle_gaps(volume)
    volume = _drop_specks(volume)

    affine = np.eye(4)
    affine[:3, 3] = lower
    return volume, affine, report


def save_phantom(path, volume, affine):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    img = nib.Nifti1Image(np.asarray(volume, np.float32), affine)
    img.header.set_xyzt_units("mm", "sec")
    img.header["descrip"] = b"right hand rule: thumb +X, index +Y, middle +Z"
    img.header.set_qform(affine, code=1)
    img.header.set_sform(affine, code=1)
    nib.save(img, str(path))
    return path


if __name__ == "__main__":
    import sys

    source = Path(sys.argv[1])
    dest = Path(sys.argv[2])
    volume, affine, report = build_phantom(source)
    save_phantom(dest, volume, affine)
    print("shape", volume.shape, "affine origin", affine[:3, 3])
    for name, item in report.items():
        if name in FINGER_NAMES.values():
            d = item["direction_mm"]
            print(f"{name:7s} direction {d[0]:7.1f} {d[1]:7.1f} {d[2]:7.1f}")
        else:
            print(name, item)
