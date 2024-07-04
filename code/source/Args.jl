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
    default = 2
    "--nmax"
    help = "Maximum number of radial basis elements"
    arg_type = Int64
    default = 2
    "--kmax"
    help = "Maximum resonance number"
    arg_type = Int64
    default = 16
    "--nbJ"
    help = "Action integral"
    arg_type = Int64
    default = 128 
    "--nbt"
    help = "Wkp integral"
    arg_type = Int64
    default = 100 
    "--nbu"
    help = "Ju computation integral"
    arg_type = Int64
    default = 100 
    "--nbv"
    help = "Jv computation integral"
    arg_type = Int64
    default = 100  
    "--nbK"
    help = "Sampling number for DF's integral"
    arg_type = Int64
    default = 100  
    "--epsLz"
    help = "Hyperparameter for semi-log evaluation of dJu/dLz"
    arg_type = Float64
    default = 0.01  
    "--errInv"
    help = "Tolerance of coordinate inversion "
    arg_type = Float64
    default = 10^(-5)  
    "--maxIter"
    help = "Max iterations for coordinate inversion "
    arg_type = Int64
    default = 50  
    "--iJmin"
    help = "Action space min index integration (for thread splitting)"
    arg_type = Int64
    default = 1
    "--iJmax"
    help = "Action space min index integration (for thread splitting)"
    arg_type = Int64
    default = 1
    "--iJmin2D"
    help = "Action space min index integration (for thread splitting) for rotation 2D integral"
    arg_type = Int64
    default = 1
    "--iJmax2D"
    help = "Action space min index integration (for thread splitting) for rotation 2D integral"
    arg_type = Int64
    default = 1
end
parsed_args = parse_args(tabargs)

 

const a = parsed_args["a"]
const c = parsed_args["c"]  
const lmax = parsed_args["lmax"]
const mmax = parsed_args["mmax"]

const nmax = parsed_args["nmax"]
const kmax = parsed_args["kmax"]
const nbJ_default = parsed_args["nbJ"]
const nbt_default = parsed_args["nbt"]
const nbu_default = parsed_args["nbu"]
const nbv_default = parsed_args["nbv"]
const nbK_default = parsed_args["nbK"]
const epsLz = parsed_args["epsLz"]
const errInverse = parsed_args["errInv"]
const maxIter_default = parsed_args["maxIter"]
const iJmin_default = parsed_args["iJmin"]
const iJmax_default = parsed_args["iJmax"]
const iJmin2d_default = parsed_args["iJmin2D"]
const iJmax2d_default = parsed_args["iJmax2D"]

@assert (a>c) "a<c is an oblate coordinate system"