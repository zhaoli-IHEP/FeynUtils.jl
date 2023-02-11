
############################################
function get_det( MM::Matrix{Basic} )::Basic
############################################

  nr, nc = size(MM)
  @assert nr == nc

  if nr == 1
    detMM = MM[1,1]
  elseif nr == 2 
    detMM = MM[1,1]*MM[2,2] - MM[1,2]*MM[2,1]
  elseif nr == 3
    detMM = MM[1,1]*MM[2,2]*MM[3,3]-MM[1,1]*MM[2,3]*MM[3,2]-MM[1,2]*MM[2,1]*MM[3,3]+
            MM[1,2]*MM[2,3]*MM[3,1]+MM[1,3]*MM[2,1]*MM[3,2]-MM[1,3]*MM[2,2]*MM[3,1]
  else
    error( "Larger nr is not expected!" )
  end # if

  return detMM

end # function get_det



####################################################
function get_adj( MM::Matrix{Basic} )::Matrix{Basic}
####################################################
  nr, nc = size( MM )
  @assert nr == nc

  adjMM = Matrix{Basic}(undef,nr,nc)
  if nr == 1

    adjMM[1,1] = Basic(1)

  elseif nr == 2

    adjMM[1,1] = MM[2,2]
    adjMM[1,2] = -MM[1,2]
    adjMM[2,1] = -MM[2,1]
    adjMM[2,2] = MM[1,1]

  elseif nr == 3

    adjMM[1,1] = -MM[2,3]*MM[3,2] + MM[2,2]*MM[3,3] 
    adjMM[1,2] = MM[1,3]*MM[3,2] - MM[1,2]*MM[3,3]
    adjMM[1,3] = -MM[1,3]*MM[2,2] + MM[1,2]*MM[2,3] 
    adjMM[2,1] = MM[2,3]*MM[3,1] - MM[2,1]*MM[3,3]
    adjMM[2,2] = -MM[1,3]*MM[3,1] + MM[1,1]*MM[3,3]
    adjMM[2,3] = MM[1,3]*MM[2,1] - MM[1,1]*MM[2,3]
    adjMM[3,1] = -MM[2,2]*MM[3,1] + MM[2,1]*MM[3,2]
    adjMM[3,2] = MM[1,2]*MM[3,1] - MM[1,1]*MM[3,2]
    adjMM[3,3] = -MM[1,2]*MM[2,1] + MM[1,1]*MM[2,2]

  else
    error( "Larger nr is not expected!" )
  end # if

  return adjMM

end # function get_adj



###############################################################
function get_dot( m1::Matrix{Basic}, m2::Matrix{Basic} )::Basic
###############################################################

  nr1, nc1 = size(m1)
  nr2, nc2 = size(m2)

  @assert nr1 == nr2 && nc1 == nc2
  nr = nr1
  nc = nc1

  if nr == 1
    vec1 = collect( view( m1, 1, : ) )
    vec2 = collect( view( m2, 1, : ) )
  elseif nc == 1
    vec1 = collect( view( m1, :, 1 ) )
    vec2 = collect( view( m2, :, 1 ) )
  else
    error( "nr, nc is not expected: "*string(nr)*", "*string(nc) )
  end # if
 
  return sum( vec1.*vec2 )

end # function get_dot











