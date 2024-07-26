
###########
## 1. SupportFun:  Partition_Generator Section 1.4.1 Generating all Partitions
## 2. SupportFun:  Partition_All
## 3. Partition_Pairs
## 4. Partition_2Perm
## 5. Partition_Indecomposable
## 6. Partition_Type_All
## 7. SupportFun: Partition_Type_eL_Location
## 8. Partition_DiagramsClosedNoLoops
## 9. Permutation_Inverse




.Partition_Pairs<- function(N){
  if (N<=2) {K2=list()
  return(K2)}
  if (N %% 2==0){
    m <- N/2
    PA<-.Partition_All(N)
    U<-PA[[1]]
    szaml<-PA[[2]]
    osz<-cumsum(szaml)
    K2<-list()
    M<-0
    for (k in 1:szaml[m]){
      K<-U[[osz[m-1]+k]]
      if (all(apply(K,1,sum)==2*rep(1,m))){
        M <- M+1
        K2[[M]]<-K
      }
    }
  }
  else {K2<-list()}
  return("K2"=K2)
}



.Partition_2Perm <-function(L){
  if (is.vector(L)){
    return("Perm_pU"=1:length(L))
  }
  n<-max(dim(L))
  h<-min(dim(L))

  ordBlock<-order(apply(L,1,sum),decreasing=T)
  pU1 <- L[ordBlock,]
  nElem<-1:n
  perm_pU<-NULL
  for (k in 1:h) {perm_pU=c(perm_pU,nElem[(pU1[k,]==1)])}
  return("Perm_pU"=perm_pU)
}


.Partition_Indecomposable <-function(L){
  if (sum(colSums(L) == rep(1,dim(L)[2]))<dim(L)[2]){
    stop("L is not a partition matrix")}
  N<-max(dim(L))
  PA<- .Partition_All(N)
  U1<-PA[[1]]
  szaml<-PA[[2]]
  IndecompK2L <- list()

  IndecompK2L[[1]]<-U1[[1]]
  M1<-NULL
  M2<-1
  osz<-cumsum(szaml)
  for (m in 2:N) {
    M<-1
    for (k in 1:szaml[m]){
      sor<-osz[m-1]
      K<-U1[[sor+k]]
      L1<- L%*%t(K)
      L2<-t(L1)%*%L1
      sL<-dim(L2)
      a1<-1
      a<-which(L2[1,]!=0)
      la<-length(a)
      la1<-length(a1)
      while (la1!=la){
        a1 <- a
        for (p in 2:length(a)){
          b<-which(L2[a[p],]!=0)
          a<-union(a,b)
        }
        la1 <- length(a1)
        la <- length(a)
      }
      if (la==sL[1]){
        M <- M+1
        M2 <- M2+1
        IndecompK2L[[M2]] <- K
      }
    }
    M1[m] <- M-1
  }

  numb_by_sizes<-c(1,M1[2:N])
  return(list("IndecompK2L"=IndecompK2L,"num_by_sizes"=numb_by_sizes))
}

#' Partitions, type and number of partitions
#'
#' Generates all partitions of \code{N} numbers and classify them by type
#'
#' @param N The (integer) number of elements to be partitioned
#' @return \code{Part.class} The list of all possible partitions given as partition matrices
#' @return \code{S_N_r}  A vector with the number of partitions of size r=1, r=2, etc. (Stirling numbers of second kind )
#' @return \code{eL_r}  A list of partition types with respect to partitions of size r=1, r=2, etc.
#' @return \code{S_r_j} Vectors of number of partitions with given types grouped by partitions of size r=1, r=2, etc.
#'
#' @examples
#' # See Example 1.18, p. 32, reference below
#' PTA<-PartitionTypeAll(4)
#' # Partitions generated
#' PTA$Part.class
#' # Partitions of size 2 includes two types
#' PTA$eL_r[[2]]
#' # Number of partitions with r=1 blocks, r=2 blocks, etc-
#' PTA$S_N_r
#' # Number of different types collected by  partitions of size r=1, r=2, etc.
#' PTA$S_r_j
#' # Partitions with size r=2, includes two types (above) each with number
#' PTA$S_r_j[[2]]
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021.  Case 1.4, p.31 and Example 1.18, p.32.
#'
#' @family Partitions
#' @export
PartitionTypeAll <- function(N){
  PA<-.Partition_All(N)
  U<-PA$U
  S_N_r<-PA$S_N_r
  sepL<-cumsum(S_N_r)

  eL_r<-list()
  S_r_j<-list()
  part_K<-list()
  part_class<-list()

  eL_r[[1]]<-c(rep(0,N-1),1)
  S_r_j[[1]]<- 1
  part_class[[1]]<-U[[1]]

  if (N==1) return(list("Part.class"=U,"S_N_r"=S_N_r,"eL_r"=eL_r,"S_r_j"=S_r_j))
  if (N==2) {
    eL_r[[N]]<-c(N,rep(0,N-1))
    S_r_j[[N]]<-1
    return(list("Part.class"=U,"S_N_r"=S_N_r,"eL_r"=eL_r,"S_r_j"=S_r_j))
  }
  for (r in 2:(N-1)){
    type_1<-matrix(0,S_N_r[r],N)

    for (k in 1:S_N_r[r]){
      type0<-matrix(0,1,N)
      K <- U[[sepL[r-1]+k]]
      sK <- apply(K,1,sum)
      C<-unique(sK)
      ia<-match(unique(sK), sK)

      if (length(ia)==1) {ind_type<-r}
      else ind_type<- diff(c(ia,r+1))

      type0[C]<- ind_type
      type_1[k,]<-type0
      part_K[[k]]<- K
    }
    y = lapply(seq_len(ncol(type_1)), function(i) type_1[,i])
    ind_type = do.call("order", c(y,list(decreasing=T)))
    type_2<-type_1[ind_type,]
    iat<-(1:nrow(type_2))[!duplicated(type_2)]
    eL_r[[r]]<-type_2[iat,]

    if (length(iat)==1) { S_r_j[[r]]<-S_N_r[[r]] }

    else S_r_j[[r]]<-diff(c(iat,S_N_r[r]+1))

    part_class<-c(part_class,part_K[ind_type])
  }
  eL_r[[N]]<-c(N,rep(0,N-1))
  S_r_j[[N]]<-1
  part_class<-c(part_class,utils::tail(U,n=1))

  return(list("Part.class"=part_class,"S_N_r"=S_N_r,"eL_r"=eL_r,"S_r_j"=S_r_j))
}





.Partition_DiagramsClosedNoLoops<- function(L){

  Kc0=list()
  N=max(dim(L));
  mi= min(dim(L));
  M=0;
  if (N%%2==0){
    m=N/2;
    U = .Partition_Indecomposable(L);
    sep=cumsum(U$num_by_sizes);

    for (k in 1:U$num_by_sizes[m]) {
      K= U$IndecompK2L[sep[m-1]+k]  ;

      if (prod(apply(K[[1]],1,sum)==rep(2,m))) {
        L1=K[[1]]%*%t(L)

        if (sum(c(L1) == 1 ) == N){
          M=M+1;
          Kc0[M]=K;
        }
      }

    }
  }
  return(Kc0)
}

############################
### 4

#' Inverse of a Permutation
#'
#' @references Gy. Terdik, Multivariate statistical methods - going beyond the linear,
#' Springer 2021, Remark 1.1, p.2
#'
#' @param permutation0 A permutation of numbers 1:n
#' @return A vector containing the inverse permutation of  \code{permutation0}
#' @family Partitions
#' @export
PermutationInv  <-  function(permutation0)
{ return(order(permutation0))}
