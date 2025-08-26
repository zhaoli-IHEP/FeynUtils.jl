


###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    is_sym_index_format( input::Basic, sym::Union{Basic, String} )::Bool

Check whether `input` is a symbol with an index like `Basic("sym1")`.
"""
function is_sym_index_format(
    input::Basic,
    sym::Union{Basic, String}
)::Bool
###################################################

  if !is_class(:Symbol,input) 
    return false
  end # if

  if !is_class(:Symbol,Basic(sym)) 
    return false
  end # if

  input_str = string(input)
  sym_str = string(sym)
  sym_len = length(sym_str)
  
  if !startswith(input_str, sym_str)
    return false
  end # if

  index_part_range  = findnext(r"[1-9]\d*", input_str, sym_len+1)
  if isnothing(index_part_range)
    return false
  end # if

  return sym_len + length(index_part_range) == length(input_str)

end # function is_sym_index_format

###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    is_loop_mom( mom::Basic )::Bool

Check whether `mom` is a loop momenta like `Basic("q1")`.
"""
is_loop_mom( mom::Basic )::Bool = is_sym_index_format(mom, "q")
###################################################

###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    is_ext_mom( mom::Basic )::Bool

Check whether `mom` is an external momenta like `Basic("K1")` or `Basic("k1")`.
Notice that `K` is for massive momentum and `k` is for massless momentum.
"""
is_ext_mom( mom::Basic )::Bool = is_sym_index_format(mom, "K") || is_sym_index_format(mom, "k")
###################################################



###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    get_sym_index( input::Basic, sym::Union{Basic, String} )::Int

Return the index for `input` respect the `sym`.
For example, `get_sym_index( Basic("x1"), "x" )` or `get_sym_index( Basic("x1"), Basic("x") )` returns `1::Int`.
Notice that `get_sym_index( Basic("x11"), "x1" )` should return `1::Int` instead of `11::Int`.
"""
function get_sym_index( 
    mom::Basic, 
    sym::Union{Basic,String} 
)::Int
###################################################

  @assert is_sym_index_format(mom, sym)
  mom_str = string(mom)
  sym_len = (length∘string)(sym)
  return parse(Int, mom_str[sym_len+1:end])

end # function get_sym_index

###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    get_loop_index( mom::Basic )::Int

Return the index of loop momenta `mom`.
For example, `get_loop_index( Basic("q1") )` returns `1::Int`.
"""
get_loop_index( mom::Basic )::Int = get_sym_index(mom, "q")
###################################################

###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    get_ext_index( mom::Basic )::Int

Return the index of external momenta `mom`.
For example, `get_loop_index( Basic("k1") )` or `get_loop_index( Basic("K1") )` returns `1::Int`.
Notice that `K` is for massive momentum and `k` is for massless momentum.
"""
get_ext_index( mom::Basic )::Int = try get_sym_index(mom, "K") catch; get_sym_index(mom, "k") end
###################################################



###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    coefficient_matrix( mom_poly_list::Vector{Basic}, mom_list::Vector{Basic} )::Matrix{Basic}

Return the coefficient matrix of momentum polynomials in `mom_poly_list` respecting to the momenta in `mom_list`.
"""
function coefficient_matrix(
    mom_poly_list::Vector{Basic}, 
    mom_list::Vector{Basic}
)::Matrix{Basic}
###################################################

  @assert all( sym -> is_class( :Symbol, sym ), mom_list )

  return reduce(
    vcat,
    (transpose ∘ map)( q -> SymEngine.coeff(mom, q, 1), mom_list )
      for mom ∈ mom_poly_list
  )

end # function coefficient_matrix



###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    get_loop_momenta( mom::Union{Basic, Vector{Basic}} )::Vector{Basic}

Find all loop momenta in `mom`.
"""
function get_loop_momenta( 
    mom::Union{Basic, Vector{Basic}}
)::Vector{Basic}
###################################################

  single_mom_list = free_symbols( mom )
  q_list = filter( is_loop_mom, single_mom_list )
  sort!( q_list; by=get_loop_index )
  return q_list

end # function get_loop_momenta



###################################################
# Created by Quan-feng Wu
# Mar 26th, 2023
"""
    get_ext_momenta( mom::Union{Basic, Vector{Basic}} )::Vector{Basic}

Find all external momenta in `mom`.
"""
function get_ext_momenta( 
    mom::Union{Basic, Vector{Basic}}
)::Vector{Basic}
###################################################

  single_mom_list = free_symbols( mom )
  k_list = filter( is_ext_mom, single_mom_list )
  sort!( k_list; by=get_ext_index )
  return k_list

end # function get_ext_momenta




