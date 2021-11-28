# Build project function
build_project () {
  target_dir=$1
  target_name=$2

  if [ "$ROOT_DIR" == "$target_dir" ]; then
    return
  fi

  cd "$target_dir" || return
  rm -rf build && mkdir build && cd build || return
  cmake .. && make || return
  cp "$target_name" "$BIN_DIR"
}

# Define and create binaries directory
BIN_DIR=$PWD/bin
if [ ! -d "$BIN_DIR" ]; then
  mkdir "$BIN_DIR"
fi

# Define root directory
ROOT_DIR=$PWD

# Build projects
build_project "$ROOT_DIR/SOR-C-reference" "sor_c"
build_project "$ROOT_DIR/SOR-MPI" "sor_mpi"
build_project "$ROOT_DIR/SOR-DSM-1" "sor_argo_1"
build_project "$ROOT_DIR/SOR-DSM-2" "sor_argo_2"
build_project "$ROOT_DIR/SOR-DSM-3" "sor_argo_3"
