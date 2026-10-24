"""Validation-only weight accounting; no nuclear model or waiting-time sampler.

Backend adapters must supply the original conditional channel probabilities
and invoke their production channel kernels. These helpers do not establish
that a backend implements the intended finite-distance transport law.
"""
from dataclasses import dataclass
import math
import random


def probabilities(values):
    values = tuple(float(v) for v in values)
    if not values or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Invalid channel probabilities')
    if not math.isclose(math.fsum(values), 1., rel_tol=0, abs_tol=1e-12):
        raise ValueError('Channel probabilities must sum to one')
    return values


def proposal(physical, importance, fraction):
    """Mix a physical branch law with a nonnegative importance guide.

    The physical component retains every possible branch. The guide only
    changes sampling; it must not supply an observable or alter a kernel.
    """
    physical = probabilities(physical)
    importance = tuple(float(v) for v in importance)
    if (len(importance) != len(physical) or
            any(not math.isfinite(v) or v < 0 for v in importance) or
            not math.isfinite(fraction) or not 0 <= fraction < 1):
        raise ValueError('Invalid importance guide or mixture fraction')
    guided = [p*g for p, g in zip(physical, importance)]
    total = math.fsum(guided)
    if not math.isfinite(total):
        raise ValueError('Importance guide overflow')
    if total == 0:
        return physical
    return probabilities([(1-fraction)*p+fraction*g/total
                          for p, g in zip(physical, guided)])


@dataclass(frozen=True)
class ChannelDraw:
    index: int
    physical_probability: float
    sampling_probability: float
    likelihood: float


def draw_channel(physical, sampling, uniform):
    """p/q applies only when holding times and all other draws are unchanged."""
    physical, sampling = probabilities(physical), probabilities(sampling)
    if (len(physical) != len(sampling) or
            any((p > 0) != (q > 0) for p, q in zip(physical, sampling)) or
            not math.isfinite(uniform) or not 0 <= uniform < 1):
        raise ValueError('Changed physical support or invalid uniform variate')
    cumulative = 0.
    selected = max(i for i, q in enumerate(sampling) if q > 0)
    for i, q in enumerate(sampling):
        cumulative += q
        if uniform < cumulative:
            selected = i
            break
    p, q = physical[selected], sampling[selected]
    return ChannelDraw(selected, p, q, p/q)


def split_weight(weight, copies):
    if (not math.isfinite(weight) or weight <= 0 or
            isinstance(copies, bool) or not isinstance(copies, int) or copies < 1):
        raise ValueError('Invalid splitting weight or replication count')
    part = weight/copies
    if part == 0:
        raise ValueError('Splitting weight underflow; never discard a branch')
    return (part,)*copies


class HistoryRewards:
    """All clones of one primary contribute to one independent observation."""
    def __init__(self):
        self._terms = {}

    def add(self, observable, weight, value):
        if (not math.isfinite(weight) or weight < 0 or
                not math.isfinite(value) or value < 0):
            raise ValueError('Invalid weighted reward')
        term = weight*value
        if not math.isfinite(term) or (weight > 0 and value > 0 and term == 0):
            raise ValueError('Weighted reward overflow or underflow')
        self._terms.setdefault(observable, []).append(term)

    def values(self):
        return {key: math.fsum(terms) for key, terms in self._terms.items()}


def competing_decay_fixture(rate, competing_rate, distance, physical, guide,
                           fraction, copies, histories, seed):
    """Analytic-method fixture only; NOT a nuclear transport backend.

    Two constant hazards compete. The first has conditional decay channels;
    the second terminates the path without a reward. Unchanged exponential
    holding times include survival to the finite boundary. Splitting occurs
    before holding-time draws, with independent continuation RNG draws.
    """
    if (not all(math.isfinite(v) and v >= 0 for v in [rate, competing_rate]) or
            not math.isfinite(distance) or distance <= 0 or
            isinstance(histories, bool) or not isinstance(histories, int) or histories < 1):
        raise ValueError('Invalid fixture parameters')
    physical = probabilities(physical)
    sampling = proposal(physical, guide, fraction)
    weights = split_weight(1., copies)
    rng = random.Random(seed)
    result = []
    total = rate+competing_rate
    if not math.isfinite(total):
        raise ValueError('Total hazard overflow')
    for _ in range(histories):
        rewards = HistoryRewards()
        for weight in weights:
            if total == 0 or rng.expovariate(total) >= distance:
                continue
            if rng.random()*total >= rate:
                continue
            draw = draw_channel(physical, sampling, rng.random())
            rewards.add(str(draw.index), weight*draw.likelihood, 1.)
        result.append(rewards.values())
    return result
