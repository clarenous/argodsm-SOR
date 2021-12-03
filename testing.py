import os

im_list = [128, 256, 512]
jm_list = [128, 256, 512]
km_list = [64, 128, 256]
node_count_list = [2, 8, 16]


def create_nx_ny_list(node_count):
    if node_count == 1:
        return [(1, 1)]
    if node_count == 2:
        return [(1, 2), (2, 1)]
    if node_count == 4:
        return [(1, 4), (2, 2), (4, 1)]
    if node_count == 8:
        return [(1, 8), (2, 4), (4, 2), (8, 1)]
    if node_count == 16:
        return [(1, 16), (2, 8), (4, 4), (8, 2), (16, 1)]
    raise Exception('unsupported ndoe_count: %d' % node_count)


def make_result_path(shape):
    result_path = "resutls/%d_%d_%d" % (shape[0], shape[1], shape[2])
    os.system('mkdir -p %s' % result_path)
    return result_path


def replace_shape(code_path, shape):
    os.system('sed -i \'s/^.*im.*$/const int im = %d;/g\' %s/sor_params.h' %
              (shape[0], code_path))
    os.system('sed -i \'s/^.*jm.*$/const int jm = %d;/g\' %s/sor_params.h' %
              (shape[1], code_path))
    os.system('sed -i \'s/^.*km.*$/const int km = %d;/g\' %s/sor_params.h' %
              (shape[2], code_path))


def replace_nx_ny(code_path, pair):
    os.system('sed -i \'s/^.*dsmNX.*$/const int dsmNX = %d;/g\' %s/sor_params.h' %
              (pair[0], code_path))
    os.system('sed -i \'s/^.*dsmNY.*$/const int dsmNY = %d;/g\' %s/sor_params.h' %
              (pair[1], code_path))


def build_target(code_path, target):
    os.system('cd %s && rm -rf build && mkdir build && cd build && cmake .. && make && mv %s ../../bin/' %
              (code_path, target))


def run_c(log_dir):
    os.system('./bin/sor_c > %s/c.out' % log_dir)


def build_and_run_c(shape):
    replace_shape("SOR-C-reference", shape)
    build_target("SOR-C-reference", "sor_c")
    run_c(make_result_path(shape))


def run_mpi(log_dir, node_count):
    os.system('mpirun -n %d ./bin/sor_mpi > %s/mpi_%d.out' %
              (node_count, log_dir, node_count))


def build_and_run_mpi(shape, node_count):
    replace_shape("SOR-MPI", shape)
    build_target("SOR-MPI", "sor_mpi")
    run_mpi(make_result_path(shape), node_count)


def run_argo(log_dir, node_count, i, pair):
    os.system('mpirun -n %d ./bin/sor_argo_%d > %s/argo_%d_%d_%d_%d.out' %
              (node_count, i, log_dir, i, node_count, pair[0], pair[1]))


def build_and_run_argo(shape, node_count):
    for i in range(1, 4):
        code_path = 'SOR-DSM-%d' % i
        replace_shape(code_path, shape)

        nx_ny_list = create_nx_ny_list(node_count)
        for pair in nx_ny_list:
            replace_nx_ny(code_path, pair)
            build_target(code_path, "sor_argo_%d" % i)
            run_argo(make_result_path(shape), node_count, i, pair)


os.system('mkdir -p bin')

for im in im_list:
    for jm in jm_list:
        for km in km_list:
            already_tested_c = False
            for node_count in node_count_list:
                print('Testing item: node_count = %d, im = %d, jm = %d, km = %d' % (node_count, im, jm, km))
                shape = (im, jm, km)
                if not already_tested_c:
                    build_and_run_c(shape)
                    already_tested_c = True
                build_and_run_mpi(shape, node_count)
                build_and_run_argo(shape, node_count)
                print()
