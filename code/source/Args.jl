using ArgParse

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--a"
    help = "Scale parameter a"
    arg_type = Float64
    default = 0.6
    "--c"
    help = "Scale parameter c"
    arg_type = Float64
    default = 0.4
    "--lmax"
    help = "Maximum harmonic number ell"
    arg_type = Int64
    default = 2
    "--mmax"
    help = "Maximum harmonic number mmax"
    arg_type = Int64
    default = 0
    "--nmax"
    help = "Maximum number of radial basis elements"
    arg_type = Int64
    default = 5
    "--kmax"
    help = "Maximum resonance number"
    arg_type = Int64
    default = 10
    # "--nbv"
    # help = "Interpolation of spherical harmonics"
    # arg_type = Int64
    # default = 1000
end
parsed_args = parse_args(tabargs)

const a = parsed_args["a"]
const c = parsed_args["c"] # Should be less than a
const lmax = parsed_args["lmax"]
const mmax = parsed_args["mmax"]

const nmax = parsed_args["nmax"]
const kmax = parsed_args["kmax"]
# const nbvInt = parsed_args["nbv"]

@assert (a>c) "a<c is an oblate coordinate system"