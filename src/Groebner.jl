

######################################
# https://math.colorado.edu/~rohi1040/expository/groebner_bases.pdf
######################################


################################
function get_monomial_weight( 
    monomial::Basic, 
    xi_list::Vector{Basic} 
)::Vector{Int64}
################################

  the_class = SymEngine.get_symengine_class(monomial)
  if the_class in [:Integer,:Rational,:Complex]
    return zeros(Int64,length(xi_list)+1)
  end # if

  @assert the_class in [:Mul,:Pow,:Symbol]
  @assert free_symbols(monomial) ⊆ xi_list

  xpt_list = map( xi -> get_exponent(monomial,xi), xi_list )

  # graded lexicographic order
  #return vcat( [sum(xpt_list)], xpt_list )
  
  # lexicographic order
  return vcat( [zero(Basic)], xpt_list )

end # function get_monomial_weight

#############################
# leading term
function get_LT( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return polynomial_expand
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  return term_list[pos]

end # function get_LT


#############################
# leading monomial e.g. power product
function get_LM( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return polynomial_expand
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  xpt_list = weight[2:end]
  return prod( xi_list .^ xpt_list )

end # function get_LM

#############################
# leading coefficient 
function get_LC( 
    polynomial::Basic, 
    xi_list::Vector{Basic} 
)::Basic
#############################

  polynomial_expand = expand(polynomial)
  the_class = SymEngine.get_symengine_class(polynomial)

  if the_class != :Add
    return subs( polynomial_expand, Dict( xi_list .=> one(Basic) ) )
  end # if

  term_list = get_add_vector_noexpand(polynomial_expand)
  weight, pos = findmax( term -> get_monomial_weight(term,xi_list), term_list )

  return subs( term_list[pos], Dict( xi_list .=> one(Basic) ) )

end # function get_LC


####################################################
# Define the function to compute the S-polynomial
function get_spoly( 
    f::Basic, 
    g::Basic, 
    xi_list::Vector{Basic} 
)::Basic
####################################################

  f_LT = get_LT( f, xi_list )
  g_LT = get_LT( g, xi_list )

  f_LT_xpt_list = get_monomial_weight( f_LT, xi_list )[2:end]
  g_LT_xpt_list = get_monomial_weight( g_LT, xi_list )[2:end]

  lcm_xpt_list = max.( f_LT_xpt_list, g_LT_xpt_list )
  lcm_fg = prod( xi_list .^ lcm_xpt_list )

  return expand( lcm_fg / f_LT * f - lcm_fg / g_LT * g )

end # function get_spoly






###########################################################
# Define the function to reduce a polynomial with respect to a Gröbner basis
function reduce_Groebner( 
    poly::Basic, 
    G_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Basic
###########################################################

  remainder = poly
  while !iszero(remainder)
    leading_term = get_LT(remainder,xi_list)
    leading_term_weight = get_monomial_weight( leading_term, xi_list )
    divided = false
    for f in G_list
      f_LT = get_LT(f,xi_list)
      f_LT_weight = get_monomial_weight( f_LT, xi_list )
      if all( f_LT_weight .<= leading_term_weight )
        factor = leading_term/f_LT
        remainder = expand( remainder - factor*f )
        divided = true
        break
      end # if
    end # for f
    if !divided
      break
    end	# if	
  end # while

  return remainder

end # function reduce_Groebner













################################################################################
# Define the function to compute the Gröbner basis using the Buchberger's algorithm
function get_Groebner_basis( 
    F_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Vector{Basic}
################################################################################

  G_list = copy(F_list)

  while true
    Gp_list = copy(G_list)
    for (p,q) in (collect∘powerset)( Gp_list, 2, 2 )
      s = get_spoly( p, q, xi_list)
      s = reduce_Groebner( s, Gp_list, xi_list )
      if !iszero(s)
        push!( G_list, s )

        G_list = setdiff( G_list, p )
        p_rem = reduce_Groebner( p, G_list, xi_list )
        if !iszero(p_rem)
          push!( G_list, p_rem )
        end # if

        G_list = setdiff( G_list, q )
        q_rem = reduce_Groebner( q, G_list, xi_list )
        if !iszero(q_rem)
          push!( G_list, q_rem )
        end # if

      end # if
    end # for one_pair

    if length(G_list) == length(Gp_list) &&
       (isempty∘setdiff)( G_list, Gp_list )
      break
    end # if
  end # while 

  while true
    Gp_list = map( x -> reduce_Groebner( x, setdiff(G_list,[x]), xi_list ), G_list )
    Gp_list = filter( !iszero, Gp_list )
    if length(G_list) == length(Gp_list) &&
       (isempty∘setdiff)( G_list, Gp_list )
      break
    else
      G_list = copy(Gp_list)
    end # if
  end # while

  for index in 1:length(G_list)
    G = G_list[index]
    G_LC = get_LC( G, xi_list )
    G_list[index] = expand(G/G_LC)
  end # for G
 
  return G_list

end # function get_Groebner_basis




######################################
# Interface to Groebner.jl
function get_Groebner_basis_v2(
    F_list::Vector{Basic}, 
    xi_list::Vector{Basic} 
)::Vector{Basic}
######################################

  xi_str_list = map( string, xi_list )
  R, new_xi_list = PolynomialRing(QQ, xi_str_list);

  polys = Vector{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}()
  for F in F_list 
    term_list = get_add_vector_expand(F)
    new_F = zero(R)
    for one_term in term_list
      the_coeff = convert( Int64, subs( one_term, Dict(xi_list .=> one(Basic)) ) )
      xpt_list = get_monomial_weight( one_term, xi_list )[2:end]
      new_F += the_coeff * prod( new_xi_list .^ xpt_list )
    end # for one_term
    push!( polys, new_F )
  end # for F

  G_list = groebner( polys, ordering=Lex() )

  new_G_list = Vector{Basic}()
  for G in G_list
    new_G = zero(Basic)
    n_term = length(G)
    for index in 1:n_term
      the_coeff = AbstractAlgebra.coeff(G,index,1) 
      one_term = term(G,index)
      xpt_list = map( x->AbstractAlgebra.degree(one_term,x), new_xi_list )

      new_term = the_coeff * prod( xi_list .^ xpt_list )
      new_G += new_term
    end # for index
    push!( new_G_list, new_G )
  end # for G

  return new_G_list

end # function get_Groebner_basis_v2

















