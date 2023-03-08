from jax import grad, jacrev, jit


def get_grad(fun, nvar):
    return jit(grad(fun, argnums=nvar))


def get_hess(fun, nvar_1, nvar_2):
    return jit(jacrev(jit(grad(fun, argnums=nvar_1)), argnums=nvar_2))
