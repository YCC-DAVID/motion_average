from sympy import Symbol, simplify
from sympy.interactive.printing import init_printing
from sympy.printing.cxx import CXX11CodePrinter
init_printing(use_unicode=False, wrap_line=False)
from sympy.matrices import Matrix, eye, zeros, ones, diag

def process_expr(ex):
    return CXX11CodePrinter().doprint(simplify(ex))

def vec_to_mat(v):
    return Matrix([[v[1], v[2], v[0]],
                   [v[4], v[5], v[3]],
                   [0, 0, 1]])

def mat_to_vec(m):
    return Matrix([m[0, 2], m[0, 0], m[0, 1], m[1, 2], m[1, 0], m[1, 1]])

if __name__ == '__main__':
    # Define Variables
    v0 = Matrix([Symbol(f'v0[{i}]') for i in range(6)])
    v1 = Matrix([Symbol(f'v1[{i}]') for i in range(6)])
    xfm = Matrix([Symbol(f'xfm[{i}]') for i in range(6)])
    P = Matrix([Symbol(f'P[{i}]') for i in range(36)])

    # Conduct dynamics
    mat_v0 = vec_to_mat(v0)
    mat_v1 = vec_to_mat(v1)

    rel01 = mat_to_vec(mat_v0.inv() @ mat_v1)
    res = rel01 - xfm

    inv0 = mat_to_vec(mat_v0.inv())
    P_mat = Matrix(P).reshape(6, 6)

    res_P = P_mat @ res
    # Generate Residual
    for i in range(6):
        print(f'residuals[{i}]={process_expr(res_P[i])};')
    print()
    print()
    ## Generate Jacobian matrix
    # dresP_dv1
    dr_dv0 = res_P.jacobian(v0)
    dr_dv1 = res_P.jacobian(v1)
    for i in range(36):
        print(f'jacobians[0][{i}]={process_expr(dr_dv0[i])};')
    print()
    for i in range(36):
        print(f'jacobians[1][{i}]={process_expr(dr_dv1[i])};')
