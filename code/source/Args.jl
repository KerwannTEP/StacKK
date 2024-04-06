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
    default = 0 # 2
    "--mmax"
    help = "Maximum harmonic number mmax"
    arg_type = Int64
    default = 0
    "--nmax"
    help = "Maximum number of radial basis elements"
    arg_type = Int64
    default = 0 # 5 
    "--kmax"
    help = "Maximum resonance number"
    arg_type = Int64
    default = 10 #16 
    "--nbJ"
    help = "Action integral"
    arg_type = Int64
    default = 10 #128
    "--nbt"
    help = "Wkp integral"
    arg_type = Int64
    default = 100 
end
parsed_args = parse_args(tabargs)

const a = parsed_args["a"]
const c = parsed_args["c"] # Should be less than a
const lmax = parsed_args["lmax"]
const mmax = parsed_args["mmax"]

const nmax = parsed_args["nmax"]
const kmax = parsed_args["kmax"]
const nbJ_default = parsed_args["nbJ"]
const nbt_default = parsed_args["nbt"]

@assert (a>c) "a<c is an oblate coordinate system"