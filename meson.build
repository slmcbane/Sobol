project('Sobol', ['cpp'])

stokhos_dep = subproject('stokhos').get_variable('stokhos_dep')

sobol_dep = declare_dependency(include_directories : include_directories('.'),
    dependencies : stokhos_dep)

if not meson.is_subproject()
    executable('gen_points', 'generate_points.cpp', dependencies : sobol_dep)
    executable('gen_random_points', 'generate_random_points.cpp', dependencies : sobol_dep)
endif
