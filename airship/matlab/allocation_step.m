function [u, achieved, gap] = allocation_step(tau, u_prev, AP)
%ALLOCATION_STEP  Share the wanted force and moment among the propellers.
%   [u, achieved, gap] = allocation_step(tau, u_prev, AP)
%   Minimises || W (A u - tau) ||^2 + 1e-8 ||u||^2 with u_lo <= u <= u_hi
%   (weights W: sideways forces and yaw first), started from the last solution u_prev.
%   achieved = A u; gap = tau - achieved (what could not be delivered).
%   Same as Allocator.solve in airship/airship/allocation.py (active-set method).
%#codegen

b = AP.A.' * (AP.W2 .* tau);
u = box_qp(AP.H, b, AP.lo, AP.hi, u_prev);
achieved = AP.A * u;
gap = tau - achieved;
end

function u = box_qp(H, b, lo, hi, u0)
% min 0.5 u'Hu - b'u  subject to lo <= u <= hi (primal active set, warm start u0)
n = numel(b);
u = min(max(u0, lo), hi);
st = zeros(n, 1);                 % 0 free, -1 at the lower bound, +1 at the upper bound
for i = 1:n
    if u(i) <= lo(i)
        st(i) = -1; u(i) = lo(i);
    elseif u(i) >= hi(i)
        st(i) = 1; u(i) = hi(i);
    end
end
for it = 1:200
    fr = st == 0;
    bd = ~fr;
    target = u;
    if any(fr)
        rhs = b(fr) - H(fr, bd) * u(bd);
        target(fr) = H(fr, fr) \ rhs;
    end
    alpha = 1.0; blk = 0; blk_st = 0;
    for i = 1:n
        if fr(i)
            d = target(i) - u(i);
            if d < 0 && target(i) < lo(i)
                a = (lo(i) - u(i)) / d;
                if a < alpha
                    alpha = a; blk = i; blk_st = -1;
                end
            elseif d > 0 && target(i) > hi(i)
                a = (hi(i) - u(i)) / d;
                if a < alpha
                    alpha = a; blk = i; blk_st = 1;
                end
            end
        end
    end
    for i = 1:n
        if fr(i)
            u(i) = u(i) + alpha * (target(i) - u(i));
        end
    end
    if blk > 0
        st(blk) = blk_st;
        if blk_st < 0
            u(blk) = lo(blk);
        else
            u(blk) = hi(blk);
        end
        continue
    end
    g = H * u - b;
    worst = 1e-12; wi = 0;
    for i = 1:n
        if st(i) == -1 && -g(i) > worst
            worst = -g(i); wi = i;
        elseif st(i) == 1 && g(i) > worst
            worst = g(i); wi = i;
        end
    end
    if wi == 0
        return
    end
    st(wi) = 0;
end
end
