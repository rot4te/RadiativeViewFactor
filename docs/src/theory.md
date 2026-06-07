# Theory

## View factor definition

The view factor ``F_{ij}`` from surface ``i`` to surface ``j`` is the fraction
of diffuse radiation leaving surface ``i`` that arrives at surface ``j``:

**3D (surface meshes):**

```math
F_{ij} = \frac{1}{A_i} \iint_{A_i} \iint_{A_j}
    \frac{\cos\theta_i \cos\theta_j}{\pi r^2} \, H_{ij} \, dA_j \, dA_i
```

**2D (curve meshes, per unit depth):**

```math
F_{ij} = \frac{1}{L_i} \int_{L_i} \int_{L_j}
    \frac{\cos\theta_i \cos\theta_j}{2r} \, H_{ij} \, dL_j \, dL_i
```

where:
- ``\theta_i`` is the angle between the outward normal at ``dA_i`` and the
  line of sight ``r_{ij}``
- ``\theta_j`` is the angle between the outward normal at ``dA_j`` and the
  reverse line of sight
- ``r`` is the distance between the two differential elements
- ``H_{ij} \in \{0, 1\}`` is the visibility function (0 = obstructed)

The factor ``2`` rather than ``\pi`` in the 2D kernel follows from integrating
the 2D radiation intensity over the hemisphere, which gives ``\pi/2`` rather
than ``\pi``.

## Reciprocity

The reciprocity relation:

```math
A_i F_{ij} = A_j F_{ji}
```

is a consequence of the symmetry of the kernel and holds exactly in the
continuous case. At the discrete (element) level, reciprocity holds to within
quadrature error; verify with [`check_reciprocity`](@ref).

## Gauss–Legendre quadrature

The isoparametric map ``\mathbf{x}(\xi, \eta)`` transforms the reference
element ``[-1,1]^2`` to physical space. The integral over one Quad8 element
becomes:

```math
\int_{A} f \, dA = \int_{-1}^{1} \int_{-1}^{1}
    f(\mathbf{x}(\xi,\eta)) \left|\frac{\partial\mathbf{x}}{\partial\xi}
    \times \frac{\partial\mathbf{x}}{\partial\eta}\right| d\xi \, d\eta
    \approx \sum_{p=1}^{n^2} w_p \, f(\mathbf{x}(\xi_p,\eta_p)) \, J_p
```

where ``J_p = |\partial_\xi \mathbf{x} \times \partial_\eta \mathbf{x}|`` is
the area Jacobian. A tensor-product ``n \times n`` Gauss–Legendre rule is used.

## Duffy transformation

For element pairs sharing a vertex at ``\mathbf{u}_0 = (u_0, v_0)`` in element
``i``'s unit square and ``\mathbf{s}_0 = (s_0, t_0)`` in element ``j``'s, the
kernel diverges as ``r \to 0``. Near the singularity:

```math
K \sim \frac{1}{r^2}, \quad r \sim \sqrt{(u-u_0)^2+(v-v_0)^2+(s-s_0)^2+(t-t_0)^2}
```

The Duffy transformation introduces a radial coordinate ``\rho`` and angular
variables ``\eta_1, \eta_2, \eta_3 \in [0,1]``. The Jacobian of the 4D
transformation is ``\rho^3``, which cancels the ``1/r^2`` divergence (since
``r \sim \rho`` near the singular corner):

```math
K \cdot dA_i \cdot dA_j \cdot |\text{Jac}| \sim \frac{1}{\rho^2} \cdot \rho^2 \cdot \rho^3 = \rho^3 \to 0
\quad \text{as } \rho \to 0
```

The transformed integrand is smooth at ``\rho = 0`` and is efficiently resolved
by standard Gauss–Legendre quadrature. The Sauter–Schwab decomposition (§5.3 of
Sauter & Schwab, 2011) provides the specific change of variables:

- **Common vertex:** 8 regions, each mapped to ``[0,1]^4``
- **Common edge:** 5 regions, each mapped to ``[0,1]^4``

## Monte Carlo estimator

The unbiased area-sampling estimator for each element pair ``(i,j)`` is:

```math
\iint K \, dA_j \, dA_i
\approx \frac{A_i \cdot A_j}{N} \sum_{k=1}^{N}
K\!\left(x_i^{(k)}, n_i^{(k)}, x_j^{(k)}, n_j^{(k)}\right) \cdot H_{ij}^{(k)}
```

where ``(x_i^{(k)}, x_j^{(k)})`` are i.i.d. uniform samples on ``A_i \times A_j``.

**Stratified sampling** subdivides the reference square into
``\lfloor\sqrt{N}\rfloor \times \lfloor\sqrt{N}\rfloor`` strata and draws one
point per stratum. For smooth integrands, stratified sampling achieves
``O(1/N)`` variance convergence rather than ``O(1/\sqrt{N})`` for plain Monte
Carlo, matching the rate of a 1D Gauss rule.

**Variance near singularities:** for the ``1/r^2`` kernel the variance of the
MC estimator is proportional to ``\iint K^2 \, dA``, which diverges at shared
edges. Infinite variance means no amount of increasing ``N`` gives reliable
convergence — use the Duffy transformation instead for such pairs.
