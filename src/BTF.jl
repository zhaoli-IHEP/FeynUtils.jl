
#-------------------------------------------------------------------------------------
# The algorithm for finding the block triangular form of the given matrix is based on
#   the reference "A graph theoretic approach to matrix inversion by partitioning"
#   the author is Frank Harary
#   publication is Numerische Mathematik 4, 128-135 (1962)
#-------------------------------------------------------------------------------------


###################################################
function generate_strong_components_list( 
    adj_mat::Matrix{Int64} 
)::Vector{Vector{Int64}}
###################################################

  n_row, n_col = size( adj_mat )
  @assert n_row == n_col
  nn = n_row

  adj_mat_plus = adj_mat + one(adj_mat)
  r_mat = map( BigInt, adj_mat_plus )^(nn-1)
  r_mat = map( entry -> iszero(entry) ? zero(Int64) : one(Int64), r_mat )

  rp_mat = transpose( r_mat )

  r_rp_mat = Matrix{Basic}( undef, nn, nn )
  for index in CartesianIndices((1:nn, 1:nn))
    r_rp_mat[index] = r_mat[index] * rp_mat[index]
  end # for index 

  strong_components_list = Vector{Int64}[]
  for row_index in 1:nn
    strong_components = sort( findall( !iszero, r_rp_mat[row_index,1:nn] ) )
    if strong_components in strong_components_list
      continue
    end # if
    push!( strong_components_list, strong_components )
  end # for row_index

  return strong_components_list

end # function generate_strong_components_list


######################################################
function generate_condensation_matrix( 
    strong_components_list::Vector{Vector{Int64}},
    star_mat::Matrix{Int64} )::Matrix{Int64}
######################################################

  n_components = length( strong_components_list )
  condensation_mat = zero( Matrix{Int64}( undef, n_components, n_components ) )
  for row_index in 1:n_components, col_index in 1:n_components
    if row_index == col_index
      continue
    end # if
    coms1 = strong_components_list[row_index]
    coms2 = strong_components_list[col_index]
    n_coms1 = length(coms1)
    n_coms2 = length(coms2)
    has_directed_path = any( index -> begin
      i1 = coms1[index[1]]
      i2 = coms2[index[2]]
      !iszero( star_mat[i1,i2] )
    end, CartesianIndices((1:n_coms1,1:n_coms2)) )
    if has_directed_path == true
      condensation_mat[row_index,col_index] = 1
    end # if
  end # for row_index, col_index

  return condensation_mat

end # function generate_condensation_matrix






#######################################################
function generate_BTF_perm( 
    the_mat::Matrix{Basic},
)::Matrix{Basic}
#######################################################

  mat_nr, mat_nc = size( the_mat )
  @assert mat_nr == mat_nc

  adj_mat = Matrix{Int64}(undef,mat_nr,mat_nc)
  for rr in 1:mat_nr, cc in 1:mat_nc
    if rr == cc || iszero(the_mat[rr,cc]) 
      adj_mat[rr,cc] = zero(Int64)
    else
      adj_mat[rr,cc] = one(Int64)
    end # if
  end # for rr, cc

  strong_components_list = generate_strong_components_list( adj_mat )

  condensation_mat = generate_condensation_matrix( strong_components_list, adj_mat )

  new_indices = Int64[]

  nr, nc = size( condensation_mat )
  tag_list = collect(1:nc)
  while !isempty( tag_list )
    nr, nc = size( condensation_mat )
    null_indices = findall( cc -> all(iszero,condensation_mat[1:nr,cc]), 1:nc )

    for one_index in null_indices
      append!( new_indices, strong_components_list[tag_list[one_index]] )
    end # for one_index
  
    condensation_mat = condensation_mat[ 1:end .∉ [null_indices], 1:end .∉ [null_indices] ]
    tag_list = tag_list[ 1:end .∉ [null_indices] ]
  end # while

  Pmat = zero( the_mat )
  for one_index in new_indices
    Pmat[new_indices[one_index],one_index] = one(Basic)
  end # for one_index

  return Pmat


end # function generate_BTF_perm











#######################################################
function generate_BTF_matrix( 
    Nmat::Matrix{Basic},
    Dmat::Matrix{Basic}
)::Tuple{ Matrix{Basic}, Matrix{Basic}, Matrix{Basic}, Vector{Int64} }
#######################################################

  mat_nr, mat_nc = size( Nmat )
  @assert mat_nr == mat_nc

  adj_mat = Matrix{Int64}(undef,mat_nr,mat_nc)
  for rr in 1:mat_nr, cc in 1:mat_nc
    if rr == cc || iszero(Nmat[rr,cc]) 
      adj_mat[rr,cc] = zero(Int64)
    else
      adj_mat[rr,cc] = one(Int64)
    end # if
  end # for rr, cc

  strong_components_list = generate_strong_components_list( adj_mat )

  condensation_mat = generate_condensation_matrix( strong_components_list, adj_mat )

  new_indices = Int64[]

  nr, nc = size( condensation_mat )
  tag_list = collect(1:nc)
  while !isempty( tag_list )
    nr, nc = size( condensation_mat )
    null_indices = findall( cc -> all(iszero,condensation_mat[1:nr,cc]), 1:nc )

    for one_index in null_indices
      append!( new_indices, strong_components_list[tag_list[one_index]] )
    end # for one_index
  
    condensation_mat = condensation_mat[ 1:end .∉ [null_indices], 1:end .∉ [null_indices] ]
    tag_list = tag_list[ 1:end .∉ [null_indices] ]
  end # while

  Pmat = zero( Nmat )
  for one_index in new_indices
    Pmat[new_indices[one_index],one_index] = one(Basic)
  end # for one_index
  Pmat⁻¹ = inv(Pmat)

  #-----------------------
  ßNmat = Pmat⁻¹*Nmat*Pmat 
  ßDmat = Pmat⁻¹*Dmat*Pmat 
  #-----------------------

  #-------------------------------
  clone_ßNmat = ßNmat
  prt_size_list = Vector{Int64}()
  while true
    nr, nc = size(clone_ßNmat) 
    @assert nr == nc
    if nr == 1
      push!( prt_size_list, 1 )
      break
    end # if

    for rr in 2:nr
      left_blk = clone_ßNmat[rr:nr,1:rr-1]
      if all( iszero, left_blk )
        nn = rr-1

        push!( prt_size_list, nn )
        clone_ßNmat = clone_ßNmat[(nn+1):end, (nn+1):end ]
        break
      end # if
    end # for rr

    if size(clone_ßNmat) == (nr, nr)
      nn = nr
      push!( prt_size_list, nn )
      clone_ßNmat = clone_ßNmat[(nn+1):end, (nn+1):end ]
      break
    end # if

  end # while
  n_prt = length(prt_size_list)

  return ßNmat, ßDmat, Pmat, prt_size_list

end # function generate_BTF_matrix 




#-------------------------------------------------------------------------------------
# Naive algorithm to construct BTF, which is used to calculate the inverse of the matrix.
#-------------------------------------------------------------------------------------

###################################
function update_indep_prt!( 
    indep_prt::Vector{Int64}, 
    flag_mat::Matrix{Int64}, 
    rem_indices::Vector{Int64}
)::Bool
###################################

  nn = (first∘size)(flag_mat)

  original_indep_prt = copy(indep_prt)

  dep_col_indices = filter( index -> any( !iszero, flag_mat[indep_prt,index] ), rem_indices )
  #dep_row_indices = filter( index -> any( !iszero, flag_mat[index,indep_prt] ), rem_indices )

  union!( indep_prt, dep_col_indices )#, dep_row_indices ) 
  sort!( indep_prt )

  return !(isempty∘setdiff)( indep_prt, original_indep_prt ) 

end # function update_indep_prt!



###################################
# find the indepdent partition, which is set of indices independent of other indices
function iterative_update_indep_prt!( 
    indep_prt::Vector{Int64}, 
    flag_mat::Matrix{Int64}, 
    rem_indices::Vector{Int64}
)::Nothing
###################################

  while true
    update_flag = update_indep_prt!( indep_prt, flag_mat, rem_indices )
    if update_flag == false
      return nothing
    end # if
  end # while

end # function iterative_update_indep_prt!



###############################
function get_block_prt_list(
    flag_mat::Matrix{Int64}, 
    rem_indices::Vector{Int64}  
)::Vector{Vector{Int64}}
###############################

  indep_prt_list = Vector{Vector{Int64}}()
  for index in rem_indices
    indep_prt = [index]
    iterative_update_indep_prt!( indep_prt, flag_mat, rem_indices )
    push!( indep_prt_list, indep_prt )
  end # for index
  unique!(indep_prt_list)
  
  min_len, min_pos = findmin( length, indep_prt_list )
  
  block_prt_list = filter( prt -> length(prt) == min_len, indep_prt_list )

  # 检验是否只有此块非零
  for block_prt in block_prt_list
    @assert all( iszero, flag_mat[block_prt,setdiff(rem_indices,block_prt)] )
  end # for index
  # 检验是否这些块之间没有重叠
  block_indices = union( block_prt_list... )
  @assert length(block_indices) == (sum∘map)( length, block_prt_list )

  # 从剩余的独立分块中再次寻找是否存在与block_prt_list中没有重叠的块
  rem_indep_prt_list = setdiff( indep_prt_list, block_prt_list ) 
  sort!( rem_indep_prt_list, by=length )
  for rem_prt in rem_indep_prt_list 
    if any( prt -> isdisjoint( prt, rem_prt ) == false, block_prt_list )
      continue
    end # if
    push!( block_prt_list, rem_prt )
  end # for indep_prt

  return block_prt_list

end # function get_block_prt_list


############################
function construct_prt(
    flag_mat::Matrix{Int64}
)::Tuple{ Vector{Int64}, Vector{Int64} }
############################

  nr, nc = size(flag_mat)
  @assert nr == nc
  nn = nr

  rem_indices = collect(1:nn)

  full_prt_list = Vector{Vector{Int64}}()
  while true
    block_prt_list = get_block_prt_list( flag_mat, rem_indices )
    append!( full_prt_list, block_prt_list )
    #@show block_prt_list
    block_indices = union( block_prt_list... )

    rem_indices = setdiff( rem_indices, block_indices )
    #@show length(rem_indices)
    if isempty(rem_indices)
      break
    end # if
  end # while

  # 反转顺序，这样返回的结果就是从上向下的上三角分块矩阵
  ordering = (reverse∘vcat)( full_prt_list... ) 
  # 因此，我们只需要分块的大小就可以了
  prt_size_list = (reverse∘map)( length, full_prt_list )

  return prt_size_list, ordering

end # function construct_prt





##########################
function block_inverse( 
    BTF_mat::Matrix{Basic}, 
    prt_size_list::Vector{Int64}
)::Matrix{Basic}
##########################

  if length(prt_size_list) == 1
    det_mat = get_det( BTF_mat )
    @assert !iszero(det_mat) 
    adj_mat = get_adj( BTF_mat )
    inv_mat = adj_mat/det_mat
    return inv_mat
  end # if

  head_prt_size = first(prt_size_list)

  A_block = BTF_mat[1:head_prt_size,1:head_prt_size]
  B_block = BTF_mat[(head_prt_size+1):end,(head_prt_size+1):end]
  C_block = BTF_mat[1:head_prt_size,(head_prt_size+1):end]

  det_A_block = get_det( A_block )
  @assert !iszero(det_A_block) 
  adj_A_block = get_adj( A_block )
  inv_A_block = adj_A_block/det_A_block

  inv_B_block = block_inverse( B_block, prt_size_list[2:end] )

  D_block = -inv_A_block*C_block*inv_B_block

  inv_mat = zero(BTF_mat)
  inv_mat[1:head_prt_size,1:head_prt_size] = inv_A_block
  inv_mat[(head_prt_size+1):end,(head_prt_size+1):end] = inv_B_block
  inv_mat[1:head_prt_size,(head_prt_size+1):end] = D_block

  #println( "obtain inverse of $(size(BTF_mat))" )
  printstyled( ".", color=:light_green ) 

  return inv_mat

end # function block_inverse






##########################
function BTF_inverse(
    the_mat::Matrix{Basic} 
)::Tuple{Matrix{Basic},Matrix{Basic}}
##########################

  nr, nc = size(the_mat)
  @assert nr == nc
  nn = nr

  flag_mat = Matrix{Int64}( undef, nn, nn )
  for rr in 1:nn, cc in 1:nn
    flag_mat[rr,cc] = iszero(the_mat[rr,cc]) ? zero(Int64) : one(Int64)
  end # for rr, cc

  cost_time = @elapsed begin
  prt_size_list, ordering = construct_prt( flag_mat )
  #@show prt_size_list ordering
  end # cost_time
  println( "[ Done construct_prt $(cost_time) sec ]" )

  cost_time = @elapsed begin
  Pmat = zeros( Basic, nn, nn )
  Pmat⁻¹ = zeros( Basic, nn, nn )
  for cc in 1:nn
    rr = ordering[cc]
    Pmat[rr,cc] = one(Basic)
    Pmat⁻¹[cc,rr] = one(Basic)
  end # for cc
  end # cost_time
  println( "[ Done $(cost_time) sec ]" )

  cost_time = @elapsed begin
  check_ordering = map( cc -> findfirst( !iszero, Pmat[:,cc] ), 1:nn )
  @assert ordering == check_ordering
  end # cost_time
  @assert Pmat*Pmat⁻¹ == one( Matrix{Basic}( undef, nn, nn ) )
  println( "[ Done $(cost_time) sec ]" )

  println( "[ Calculate BTF_mat ]" )
  flush(stdout)
  cost_time = @elapsed begin
  BTF_mat = Pmat⁻¹*the_mat*Pmat
  end # cost_time
  println( "[ Done multiplication $(cost_time) sec ]" )

  cost_time = @elapsed begin
  check_BTF_mat = the_mat[ordering,ordering]
  check_mat = rational_function_simplify( BTF_mat - check_BTF_mat )
  @assert all( iszero, check_mat )

  for index in 1:length(prt_size_list)
    prt_size = prt_size_list[index]
    header = sum( prt_size_list[1:(index-1)] ) 
    check = all( iszero, BTF_mat[(header+prt_size+1):end,(header+1):(header+prt_size)] )
    #@show index check
    @assert check == true
  end # for prt_size

  end # cost_time
  println( "[ Done check $(cost_time) sec ]" )


  println( "[ Calculate block inverse ]" )
  #@show prt_size_list
  cost_time = @elapsed begin
  inv_BTF_mat = block_inverse( BTF_mat, prt_size_list )
  end # cost_time
  println( "[ Done $(cost_time) sec ]" )

  inv_ordering = map( cc -> findfirst( !iszero, Pmat⁻¹[:,cc] ), 1:nn )
  #@show inv_ordering
  inv_the_mat = inv_BTF_mat[inv_ordering,inv_ordering]
  
  println( "Check the_mat*inv_the_mat" )
  cost_time = @elapsed begin
  check_mat = rational_function_simplify( the_mat*inv_the_mat - one(the_mat) )
  end # cost_time
  println( "[ Done $(cost_time) sec ]" )
  flush(stdout)
  @assert all( iszero, check_mat )

  println( "[ Simplify the inverse matrix ]" )
  cost_time = @elapsed begin
  num_inv_the_mat, den_inv_the_mat = numer_denom_simplify( inv_the_mat )
  end # cost_time
  println( "[ Done $(cost_time) sec ]" )
  flush(stdout)

  return num_inv_the_mat, den_inv_the_mat

end # function BTF_inverse







