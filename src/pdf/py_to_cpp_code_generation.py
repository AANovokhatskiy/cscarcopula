import sympy as sp

from sympy.printing import ccode
from sympy import symbols, cse, numbered_symbols, Symbol
import re


def generator(copula_name):
    t, r = sp.symbols('t r')
    if copula_name == 'Gumbel':
        sp_generator = (-sp.log(t))**r
        sp_inverse_generator = sp.exp(-t**(1/r))
        transform = sp.Min(r**2 + 1, 40)
    else:
        raise ValueError(f'code generation is not implemented for {copula_name} copula')
    
    return sp_generator, sp_inverse_generator, transform
 

def sp_pdf(dim, sp_generator, sp_inverse_generator):
    u = sp.symbols('u0:%d'%(dim))
    t = sp.symbols('t')
    params = [sp_generator.subs([(t, x)]) for x in u]
    diff_inverse_generator = sp.together(sp_inverse_generator.diff((t, dim)))
    diff_generator = sp.together(sp_generator.diff(t, 1))
    func = diff_inverse_generator.subs([(t, sum(params))])
    for x in u:
        func = func * diff_generator.subs([(t, x)])
    return sp.powsimp(func)


# def sympy_to_c(function_name, sympy_function):
#     tmpsyms = numbered_symbols("tmp")
#     symbols, simple = cse(sympy_function, symbols=tmpsyms)
#     symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
#     varstring=", ".join("double "+ x for x in symbolslist)
#     c_code = "double "+str(function_name)+" ("+varstring+")\n"
#     c_code +=  "{\n"
#     for s in symbols:
#         #print s
#         c_code +=  "  double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
#     c_code +=  "  result = " + ccode(simple[0])+";\n"
#     c_code +=  "  return result;\n"
#     c_code += "}\n"
#     return c_code

def sympy_to_c_pdf_double_r_vector_vector_u(function_name, sympy_function):
    tmpsyms = numbered_symbols("tmp")

    symbols, simple = cse(sympy_function, symbols=tmpsyms)
    symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
    c_code = f"std::vector<double> {function_name}(const std::vector<std::vector<double>>& u, const double r)\n"
    c_code +=  "{\n"
    c_code += "    size_t n = u.size();\n"
    c_code += "    std::vector<double> result(n);\n \n"
    c_code += "    for (size_t i = 0; i < n; ++i) {\n"
    for j in range(0, len(symbolslist) - 1):
        c_code += f"        double u{j} = u[i][{j}];\n"
    for s in symbols:
        c_code +=  "        double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
    c_code +=  "        result[i] = " + ccode(simple[0])+";\n"
    c_code +=  "    }\n"
    c_code +=  "  return result;\n"
    c_code += "}\n"
    return c_code

def sympy_to_c_pdf_double_r_vector_u(function_name, sympy_function):
    tmpsyms = numbered_symbols("tmp")

    symbols, simple = cse(sympy_function, symbols=tmpsyms)
    symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
    c_code = f"double {function_name}(const std::vector<double>& u, const double r)\n"
    c_code +=  "{\n"
    c_code += "    double result = 0.0;\n \n"
    for j in range(0, len(symbolslist) - 1):
        c_code += f"    double u{j} = u[{j}];\n"
    for s in symbols:
        c_code +=  "    double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
    c_code +=  "        result += " + ccode(simple[0])+";\n"
    c_code +=  "  return result;\n"
    c_code += "}\n"
    return c_code

def sympy_to_c_pdf_vector_r_vector_vector_u(function_name, sympy_function):
    tmpsyms = numbered_symbols("tmp")

    r = sp.symbols('r')
    aux_r_name = "aux_r_expr"
    sympy_function = sympy_function.subs(r, aux_r_name)

    symbols, simple = cse(sympy_function, symbols=tmpsyms)
    symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
    c_code = f"std::vector<double> {function_name}(const std::vector<std::vector<double>>& u, const std::vector<double>& r)\n"
    c_code +=  "{\n"
    c_code += "    size_t n = u.size();\n"
    c_code += "    std::vector<double> result(n);\n \n"
    c_code += "    for (size_t i = 0; i < n; ++i) {\n"
    for j in range(0, len(symbolslist) - 1):
        c_code += f"        double u{j} = u[i][{j}];\n"
    for s in symbols:
        s0 = re.sub(aux_r_name, "r[i]", ccode(s[0]))
        s1 = re.sub(aux_r_name, "r[i]", ccode(s[1]))
        c_code +=  "        double " + s0 + " = " + s1 + ";\n"
    s2 = re.sub(aux_r_name, "r[i]", ccode(simple[0]))
    c_code +=  "        result[i] = " + s2 +";\n"
    c_code +=  "    }\n"
    c_code +=  "  return result;\n"
    c_code += "}\n"
    return c_code

def sympy_to_c_double_r(function_name, sympy_function):
    tmpsyms = numbered_symbols("tmp")
    symbols, simple = cse(sympy_function, symbols=tmpsyms)
    symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
    c_code = f"double {function_name}(const double r)\n"
    c_code +=  "{\n"
    c_code += "    double result = 0.0;\n"
    for s in symbols:
        c_code +=  "  double " + ccode(s[0]) + " = " + ccode(s[1]) + ";\n"
    c_code +=  "  result = " + ccode(simple[0])+";\n"
    c_code +=  "  return result;\n"
    c_code += "}\n"
    return c_code

def sympy_to_c_vector_r(function_name, sympy_function):
    tmpsyms = numbered_symbols("tmp")
    r = sp.symbols('r')
    aux_r_name = "aux_r_expr"
    sympy_function = sympy_function.subs(r, aux_r_name)

    symbols, simple = cse(sympy_function, symbols=tmpsyms)
    symbolslist = sorted(map(lambda x: str(x), list(sympy_function.atoms(Symbol))))
    c_code = f"std::vector<double> {function_name}(const std::vector<double>& r)\n"
    c_code +=  "{\n"
    c_code += "    size_t n = r.size();\n"
    c_code += "    std::vector<double> result(n);\n \n"
    c_code += "    for (size_t i = 0; i < n; ++i) {\n"
    for s in symbols:
        s0 = re.sub(aux_r_name, "r[i]", ccode(s[0]))
        s1 = re.sub(aux_r_name, "r[i]", ccode(s[1]))
        c_code +=  "        double " + s0 + " = " + s1 + ";\n"
    s2 = re.sub(aux_r_name, "r[i]", ccode(simple[0]))
    c_code +=  "        result[i] = " + s2 +";\n"
    c_code +=  "    }\n"
    c_code +=  "  return result;\n"
    c_code += "}\n"
    return c_code

def h_header(copula_name, dim):
    c_code = f"//{copula_name}\n"
    c_code += f"//{dim}\n \n"
    c_code += "#pragma once\n \n"
    c_code += "#include <vector>\n"
    c_code += "#include <cmath>\n \n"
    c_code += "std::vector<double> pdf(const std::vector<std::vector<double>>& u, const double r);\n \n"
    c_code += "std::vector<double> pdf(const std::vector<std::vector<double>>& u, const std::vector<double>& r);\n \n"
    c_code += "double pdf(const std::vector<double>& u, const double r);\n \n"
    c_code += "double transform(const double r);\n \n"
    c_code += "std::vector<double> transform(const std::vector<double>& r);\n \n"
    return c_code



if __name__ == '__main__':
    dim = 9
    copula_name = 'Gumbel'
    gen, inv_gen, transform = generator('Gumbel')
    pdf = sp_pdf(dim, gen, inv_gen)
    header = h_header(copula_name, dim)
    # cpp_fun = sympy_to_c('pdf', pdf )
    cpp_fun_double_r_vector_vector_u = sympy_to_c_pdf_double_r_vector_vector_u('pdf', pdf )
    cpp_fun_vector_r_vector_vector_u = sympy_to_c_pdf_vector_r_vector_vector_u('pdf', pdf )
    cpp_fun_vector_r_vector_u = sympy_to_c_pdf_double_r_vector_u('pdf', pdf )
    cpp_fun_double_r = sympy_to_c_double_r('transform', transform)
    cpp_fun_vector_r = sympy_to_c_vector_r('transform', transform)

    with open("generated_pdf.h", "w") as output:
        output.write(header)
    with open("generated_pdf.cpp", "w") as output:
        output.write("""#include "generated_pdf.h" """)
        output.write("\n \n")
        # output.write(cpp_fun)
        # output.write("\n")
        output.write(cpp_fun_double_r_vector_vector_u)
        output.write("\n")
        output.write(cpp_fun_vector_r_vector_vector_u)
        output.write("\n")
        output.write(cpp_fun_vector_r_vector_u)
        output.write("\n")
        output.write(cpp_fun_double_r)
        output.write("\n")
        output.write(cpp_fun_vector_r)
