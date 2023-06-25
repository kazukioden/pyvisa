program wile
!wile.f90に磁場を追加しただけ
use EVLCG_INT
use WRCRN_INT
use WRRRN_INT
	implicit none
integer  i,j,n,p,q,m,nm,mh,nso,time,ITMAX,alert,k,lp,lps,tansaku,zu
complex(8),dimension(4) :: cnumber,v1
real(8) A_d(4), B_d(4),IP(4),sig(8),sha(2),thetax,tpop(4,3,10),ypo(4,100),potan(4,3),tan(4,3,100),ych2(100),poselect(4,3),youb,pomod(4,3)
real(8) P_so,y,z,phi1,phi2,phi3,EVALDRE(8),ER(4),EI(4),po(4,3),firstpo(4,3),poswap(3),ftol,yp(4),firstyp(4),ych(100),ype(4),ype2(4)
real(8),allocatable :: pops(:,:,:),ypos(:,:)
real(8),parameter :: pi =3.14159265359
complex(8) check
complex(8),dimension(4) :: l,a,b,ad,bd
complex(8),dimension(8,8) :: c
real(8),dimension(4) ::ar,ai,br,bi,clr,cli,phi,theta
real(8),dimension(4,4) :: d_,e_
real(8),dimension(8,8) :: x,xt !tはtildeの略
complex(8),dimension(8,8) :: ct, cd !dはdashの略
complex(8),dimension(8,8) :: cdt
complex(8),dimension(8,8) :: U,Udag !dagはdaggerの略
complex(8),dimension(8,8) :: lambda,test,Sd,S
complex(8),dimension(8) :: EVAL,EVALD
complex(8),dimension(8,8) :: EVEC,SL,SLd,expmiphi,expiphi,tau
!毎回違う乱数生成のための記述
integer :: seedsize
integer,allocatable :: seed(:)
seedsize=1000000
call random_seed(size=seedsize)
allocate(seed(seedsize))
do i=1,seedsize
 call system_clock(count=seed(i))
end do
 seed(:)= 1378046362
call random_seed(put=seed(:))
 time=0
 alert=0

open(3,file='tokuitejiba.csv',status='replace')
open(5,file='tokuitejiba1e.csv',status='replace')
open(7,file='tokuitejiba112.csv',status='replace')
open(33,file='tokuitejiba123.csv',status='replace')
open(35,file='tokuitejiba113.csv',status='replace')

call random_number(d_)
d_=chirasu2(d_)
call random_number(e_)
e_=chirasu2(e_)
call random_number(ar)
ar=chirasu1(ar)
call random_number(ai)
ai=chirasu1(ai)
call random_number(br)
br=chirasu1(br)
call random_number(bi)
bi=chirasu1(bi)
call random_number(clr)
clr=chirasu1(clr)
call random_number(cli)
cli=chirasu1(cli)

!ここではPsoと磁場を調整　ここではzuという整数が18まではpsoが増え、その後は磁場が増えるようになってる　あとはwile.f90と同じと言いたいところだけど、338行目少し変えたわ　一つ前の条件のワイル特異点からエネルギーのズレが最小になるようになるものを選んでる
nso=0
mh=1

do zu=1,2909!2946

if(zu<18) then
nso=nso+1
else if(zu>18) then
mh=mh+1
end if


!do nso=1,176
P_so=0.01d0*(-1.0d0+dble(nso))


thetax=dble(mh-1)*pi/20000.0d0



!ランダム生成






do i=1,4
a(i)=cmplx(ar(i),ai(i))  
b(i)=cmplx(br(i),bi(i))  
end do
ad=a
bd=b
do i=1,4
v1(i)=sqrt(dconjg(ad(i))*ad(i)+dconjg(bd(i))*bd(i))
if (abs(v1(i))<10**(-3)) then
a(i)=0.0d0
b(i)=0.0d0
else
a(i)=ad(i)/v1(i)
b(i)=bd(i)/v1(i)
end if
enddo

do i=1,4
 do j=1,4
 c(2*j-1,2*i-1)=cmplx(d_(j,i),e_(j,i))*a(i)
 c(2*j,2*i-1)=cmplx(d_(j,i),e_(j,i))*b(i)
 c(2*j-1,2*i)=-dconjg(c(2*j-1,2*i-1))
 c(2*j,2*i)=dconjg(c(2*j,2*i-1))
 end do
end do
print *,"size=",size(c,2)
!cのシュミットの直交化
ct=gs(c,8)
!print *,ct
sha=shape(c)
!print *,'shapec=',sha(2)


!チェック
!do i=1,8
 !do j=1,8
   !check=dot_product1(dconjg(ct(1:8,j)),ct(1:8,i))
  !print *,"check1=",check
 !enddo
!enddo
!print *,"check1=",check
!!!!!('ct',ct,8,8,8)
!実部のシュミットの直交化
x=dreal(ct)
do i=1,4
 do j=1,4
  x(2*j-1,2*i-1)=x(2*j-1,2*i-1)
  x(2*j,2*i-1)=0
  x(2*j-1,2*i)=0
  x(2*j,2*i)=x(2*j-1,2*i-1)
  end do
end do
xt=gsr(x,8)
!call WRRRN('xt',xt,8,8,8)
!チェック
!check=0
!do i=1,8
 !do j=1,8
  
   !check=dot_product2((xt(1:8,j)),xt(1:8,i))
  
 ! print *,"check2=",check
 !enddo
!enddo
!print *,"check2=",check

!P_soパラメータを含んだcd(j,i)
do i=1,8
do j=1,8
cd(j,i)=xt(j,i)+P_so*(ct(j,i)-xt(j,i))
end do 
end do
!cd(j,i)の直交化
cdt=gs(cd,8)
!print *,cdt
!チェック
!do i=1,8
 !do j=1,8
  ! check=dot_product1(dconjg(cdt(1:8,j)),cdt(1:8,i))
  !print *,"check3=",check
 !enddo
!enddo
!do i=1,8
 !do j=1,8
  ! check=dot_product1(cdt(1:8,j),dconjg(cdt(1:8,i)))
  !print *,"check4=",check
 !enddo
!enddo
!ユニタリ行列Uを構成
U=cdt
do i=1,4
 do j=1,4
 U(2*j-1,2*i-1)=cdt(2*j-1,2*i-1)
 U(2*j,2*i-1)=cdt(2*j,2*i-1)
 U(2*j-1,2*i)=-dconjg(U(2*j,2*i-1))
 U(2*j,2*i)=dconjg(U(2*j-1,2*i-1))
 end do
end do
!!!!!('U',U,8,8,8)

!Uの複素転置行列を求める
Udag=dconjg(transpose(U))
!Uのユニタリー性をテスト

!!!!!('test1',test,8,8,8)



!lambdaを構成


do i=1,4
cnumber(i)=cmplx(clr(i),cli(i))
v1(i)=sqrt(dconjg(cnumber(i))*cnumber(i))
cnumber(i)=cnumber(i)/v1(i)
enddo

do i=1,4
lambda(2*i-1,2*i-1)=cnumber(i)
lambda(2*i,2*i)=cnumber(i)
end do
!!!!!('lambda',lambda,8,8,8)



!Sを構成
Sd=matmul(lambda,dconjg(transpose(U)))
S=matmul(U,Sd)
!!!!!('S',S,8,8,8)
!print *,S(1,2)
!Sのユニタリー性をテスト
!test=matmul(dconjg(transpose(S)),S)
!!!!!('test2',test,8,8,8)
!S*Sを計算してみる
!test=matmul(dconjg(S),S)
!!!!!('S*S',test,8,8,8)



if (time==0) then
call random_number(po)
po=chirasu2(po)
po=po*pi
do i=1,4
yp(i)=andreev(po(i,:),S,thetax)
end do
ftol=0.0000001d0
ITMAX=200000

call  amoeba(po,yp,4,3,3,ftol,andreev,S,ITMAX,thetax)
firstpo=po
firstyp=yp
 do i=1,3
  do j=i+1,4
   if(dabs(yp(i))<=dabs(yp(j))) then
   else
    y=yp(i)
    yp(i)=yp(j)
    yp(j)=y
    poswap=po(i,:)
    po(i,:)=po(j,:)
    po(j,:)=poswap
   end if
  end do
 end do

if(yp(1)>0.001d0) then
call WRRRN('firstpo',firstpo,4,3,4)
 call WRRRN('firstyp',firstyp,4,1,4)
 print *,"seed=",seed,"P_so=",P_so,"time=",time
 stop "hazure"
end if
 endif


!if(yp(1)>0.01d0) then
!alert=alert+1
!endif

!if(alert>10) then
!call WRRRN('firstpo',firstpo,4,3,4)
! call WRRRN('firstyp',firstyp,4,1,4)
! print *,"seed=",seed,"P_so=",P_so,"time=",time
! stop "kowareta"
!end if


 !if(yp(1)<0.01d0.and.time>5) then
 !print *,"seed=",seed,"P_so=",P_so
 !stop "kowareta"
 !endif
 if(time/=0) then
 ITMAX=20000
 ftol=0.0000001d0
 

do tansaku=1,100
 call random_number(tan(:,:,tansaku))
 tan(:,:,tansaku)=chirasu2(tan(:,:,tansaku))
 potan=po+pi*0.01d0*tan(:,:,tansaku)








do i=1,4
yp(i)=andreev2(potan(i,:),S,thetax)
end do

call  amoeba2(potan(:,:),yp(:),4,3,3,ftol,andreev2,S,ITMAX,thetax)
do i=1,3
  do j=i+1,4
   if(dabs(yp(i))>dabs(yp(j))) then
    y=yp(i)
    yp(i)=yp(j)
    yp(j)=y
    poswap=potan(i,:)
    potan(i,:)=potan(j,:)
    potan(j,:)=poswap
   end if
  end do
 end do

 do i=1,4
 ype(i)=andreev(potan(i,:),S,thetax)
 end do

 do i=1,4
   do j=i+1,4
  if(dabs(ype(i))>dabs(ype(j))) then
    y=ype(i)
    ype(i)=ype(j)
    ype(j)=y
  end if
   end do
 end do
    ych(tansaku)=ype(1)

   


 


 
 if (tansaku/=1 .and.dabs(ych(tansaku)-youb)<dabs(ych(1)-youb)) then
 ych(1)=ych(tansaku)
 poselect=potan
 end if
end do
po=poselect

 do i=1,4
 ype2(i)=andreev2(poselect(i,:),S,thetax)
 end do

end if
!if(ype2(1)>0.01d0) then
!stop 'kowaretapsode'
!end if
pomod=po
 do i=1,3
 if (pomod(1,i)>pi) then
  pomod(1,i)=pomod(1,i)-2*pi
 else if (pomod(1,i)<-pi) then
  pomod(1,i)=pomod(1,i)+2*pi
 end if
 end do
write(3,*)thetax,",",P_so,",",ych(1),",",ype2(1),",",po(1,1),",",po(1,2),",",po(1,3)
write(5,*)thetax/pi,",",ych(1)!,",",ype2(1),",",po(1,1),",",po(1,2),",",po(1,3)
write(7,*)po(1,1)/pi,",",po(1,2)/pi!,",",po(1,3)
write(33,*)po(1,2)/pi,",",po(1,3)/pi!,",",po(1,3)
write(35,*)po(1,1)/pi,",",po(1,3)/pi!,",",po(1,3)
time=time+1
youb=ych(1)
end do


call random_number(tpop)
tpop=chirasu3(tpop)*pi 

 
 lp=0


do k=1,10
 do i=1,4
  ypo(i,k)=andreev(tpop(i,:,k),S,thetax)
 end do
call  amoeba(tpop(:,:,k),ypo(:,k),4,3,3,ftol,andreev,S,ITMAX,thetax) 
end do

do k=1,10
 do i=1,3
  do j=i+1,4
   if(dabs(ypo(i,k))<=dabs(ypo(j,k))) then
   else
    y=ypo(i,k)
    ypo(i,k)=ypo(j,k)
    ypo(j,k)=y
    poswap=tpop(i,:,k)
    tpop(i,:,k)=tpop(j,:,k)
    tpop(j,:,k)=poswap
   end if
  end do
 end do
end do


do k=1,10
if (dabs(ypo(1,k))<0.001d0) then
lp=lp+1
endif
enddo
allocate(ypos(4,lp))
allocate(pops(4,3,lp))

lps=0
do k=1,10
if (dabs(ypo(1,k))<0.001d0) then
lps=lps+1
ypos(:,lps)=ypo(:,k)
pops(:,:,lps)=tpop(:,:,k)
end if
end do




call WRRRN('po',po,4,3,4)




call WRRRN('firstpo',firstpo,4,3,4)
 call WRRRN('firstyp',firstyp,4,1,4)
 print *,"seed=",seed,"P_so=",P_so,"time=",time

do k=1,lp
call WRRRN('pops(:,:,k)',pops(:,:,k),4,3,4)
call WRRRN('ypos(:,k)',ypos(:,k),4,1,4)
end do


 contains
   SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,andreev,Stot,ITMAX,thetax)
   IMPLICIT NONE
   INTEGER  iter,mp,ndim,np,k,ITMAX
   real(8)  ftol,p(mp,np),y(mp),andreev,thetax
   complex(8),dimension(8,8) ::  Stot
   integer,PARAMETER :: NMAX=3
   real(8),parameter :: TINY=1.0d-9

   INTEGER :: i,ihi,ilo,inhi,j,m,n
 double precision :: rtol,sum,swap,ysave,ytry,psum(NMAX)!,amotry
   iter=0
    do n=1,ndim 
     sum=0.0d0 
      do m=1,ndim+1
       sum=sum+p(m,n)
      end do
     psum(n)=sum
    end do

    do k=1,ITMAX
     ilo=1
      if (y(1).gt.y(2)) then
       ihi=1
       inhi=2
      else
       ihi=2
       inhi=1
      endif
     do i=1,ndim+1
      if(y(i).le.y(ilo)) ilo=i
      if(y(i).gt.y(ihi)) then
       inhi=ihi
       ihi=i
      else if(y(i).gt.y(inhi)) then
       if(i.ne.ihi) inhi=i
      endif
     end do
     rtol=2.0d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo))+TINY)
      if (rtol.lt.ftol) then
       swap=y(1)
       y(1)=y(ilo)
       y(ilo)=swap
        do n=1,ndim
         swap=p(1,n)
         p(1,n)=p(ilo,n)
         p(ilo,n)=swap
        end do
       return
      endif
     !if (iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
     iter=iter+2
     ytry=amotry(p,y,psum,mp,np,ndim,andreev,Stot,ihi,-1.0d0,thetax)
     if (ytry.le.y(ilo)) then

      ytry=amotry(p,y,psum,mp,np,ndim,andreev,Stot,ihi,2.0d0,thetax)
     else if (ytry.ge.y(inhi)) then

      ysave=y(ihi)
      ytry=amotry(p,y,psum,mp,np,ndim,andreev,Stot,ihi,0.5d0,thetax)
      if (ytry.ge.ysave) then
       do i=1,ndim+1
        if(i.ne.ilo)then
         do j=1,ndim
          psum(j)=0.5d0*(p(i,j)+p(ilo,j))
          p(i,j)=psum(j)
         end do
        y(i)=andreev(psum,Stot,thetax)
        endif
       end do
      iter=iter+ndim
       do n=1,ndim 
        sum=0.0d0 
         do m=1,ndim+1
          sum=sum+p(m,n)
         end do
        psum(n)=sum
       end do
      endif
     else
       iter=iter-1
     endif
    end do !k
   END SUBROUTINE amoeba


   FUNCTION amotry(p,y,psum,mp,np,ndim,andreev,Stot,ihi,fac,thetax)
   IMPLICIT NONE
   INTEGER :: ihi,mp,ndim,np
   real(8)  amotry,fac,p(mp,np),psum(np),y(mp),andreev,thetax
   complex(8)  Stot(8,8)
   integer,PARAMETER :: NMAX=3

   INTEGER :: j
   real(8)  fac1,fac2,ytry,ptry(NMAX)
   fac1=(1.0d0-fac)/dble(ndim)
   fac2=fac1-fac
    do j=1,ndim
     ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
    end do
   ytry=andreev(ptry,Stot,thetax)
    if (ytry.lt.y(ihi)) then
     y(ihi)=ytry
      do j=1,ndim
       psum(j)=psum(j)-p(ihi,j)+ptry(j)
       p(ihi,j)=ptry(j)
      end do
    endif
   amotry=ytry
   return
   END FUNCTION amotry


   SUBROUTINE amoeba2(p,y,mp,np,ndim,ftol,andreev2,Stot,ITMAX,thetax)
   IMPLICIT NONE
   INTEGER  iter,mp,ndim,np,k,ITMAX
   real(8)  ftol,p(mp,np),y(mp),andreev2,thetax
   complex(8),dimension(8,8) ::  Stot
   integer,PARAMETER :: NMAX=3
   real(8),parameter :: TINY=1.0d-9

   INTEGER :: i,ihi,ilo,inhi,j,m,n
 real(8)  rtol,sum,swap,ysave,ytry,psum(NMAX)!,amotry2
   iter=0
    do n=1,ndim 
     sum=0.0d0 
      do m=1,ndim+1
       sum=sum+p(m,n)
      end do
     psum(n)=sum
    end do

    do k=1,ITMAX
     ilo=1
      if (y(1).gt.y(2)) then
       ihi=1
       inhi=2
      else
       ihi=2
       inhi=1
      endif
     do i=1,ndim+1
      if(y(i).le.y(ilo)) ilo=i
      if(y(i).gt.y(ihi)) then
       inhi=ihi
       ihi=i
      else if(y(i).gt.y(inhi)) then
       if(i.ne.ihi) inhi=i
      endif
     end do
     rtol=2.0d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo))+TINY)
      if (rtol.lt.ftol) then
       swap=y(1)
       y(1)=y(ilo)
       y(ilo)=swap
        do n=1,ndim
         swap=p(1,n)
         p(1,n)=p(ilo,n)
         p(ilo,n)=swap
        end do
       return
      endif
     !if (iter.ge.ITMAX) pause 'ITMAX exceeded in amoeba'
     iter=iter+2
     ytry=amotry2(p,y,psum,mp,np,ndim,andreev2,Stot,ihi,-1.0d0,thetax)
     if (ytry.le.y(ilo)) then

      ytry=amotry2(p,y,psum,mp,np,ndim,andreev2,Stot,ihi,2.0d0,thetax)
     else if (ytry.ge.y(inhi)) then

      ysave=y(ihi)
      ytry=amotry2(p,y,psum,mp,np,ndim,andreev2,Stot,ihi,0.5d0,thetax)
      if (ytry.ge.ysave) then
       do i=1,ndim+1
        if(i.ne.ilo)then
         do j=1,ndim
          psum(j)=0.5d0*(p(i,j)+p(ilo,j))
          p(i,j)=psum(j)
         end do
        y(i)=andreev2(psum,Stot,thetax)
        endif
       end do
      iter=iter+ndim
       do n=1,ndim 
        sum=0.0d0 
         do m=1,ndim+1
          sum=sum+p(m,n)
         end do
        psum(n)=sum
       end do
      endif
     else
       iter=iter-1
     endif
    end do !k
   END SUBROUTINE amoeba2


   FUNCTION amotry2(p,y,psum,mp,np,ndim,andreev2,Stot,ihi,fac,thetax)
   IMPLICIT NONE
   INTEGER :: ihi,mp,ndim,np
   real(8)  amotry2,fac,p(mp,np),psum(np),y(mp),andreev2,thetax
   complex(8)  Stot(8,8)
   integer,PARAMETER :: NMAX=3

   INTEGER :: j
   real(8)  fac1,fac2,ytry,ptry(NMAX)
   fac1=(1.0d0-fac)/dble(ndim)
   fac2=fac1-fac
    do j=1,ndim
     ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
    end do
   ytry=andreev2(ptry,Stot,thetax)
    if (ytry.lt.y(ihi)) then
     y(ihi)=ytry
      do j=1,ndim
       psum(j)=psum(j)-p(ihi,j)+ptry(j)
       p(ihi,j)=ptry(j)
      end do
    endif
   amotry2=ytry
   return
   END FUNCTION amotry2
function andreev(phid,S,thetax) result(wild)
implicit none
complex(8),intent(in) :: S(8,8)
real(8),intent(in) :: phid(3),thetax
real(8) wild,kikaku(8)
integer  i,j,n,p,q,m,nm,mh
complex(8),dimension(4) :: cnumber,v1
real(8) A_d(4), B_d(4),IP(4),sig(8),sha(2)
real(8) P_so,y,z,phi1,phi2,phi3,EVALDRE(8),ER(4),EI(4)
real(8),parameter :: pi =3.14159265359
complex(8) check
complex(8),dimension(4) :: l,a,b,ad,bd
complex(8),dimension(8,8) :: c
real(8),dimension(4) ::ar,ai,br,bi,clr,cli,phi,theta
real(8),dimension(4,4) :: d_,e_
real(8),dimension(8,8) :: x,xt !tはtildeの略
complex(8),dimension(8,8) :: ct, cd !dはdashの略
complex(8),dimension(8,8) :: cdt
complex(8),dimension(8,8) :: U,Udag !dagはdaggerの略
complex(8),dimension(8,8) :: lambda,test,Sd
complex(8),dimension(8) :: EVAL,EVALD,EDA,EVALa
complex(8),dimension(8,8) :: EVEC,SL,SLd,expmiphi,expiphi,tau

!SLを構成
!tauを構成し、かける
!call WRCRN('S',S,8,8,8)
!test=matmul(dconjg(transpose(S)),S)
!call WRCRN('test',test,8,8,8)
do i=1,4
theta(i)=thetax
end do
!do i=1,4
 !tau(2*i-1,2*i-1)=dexp(cmplx(0.0d0,theta(i)))
  !tau(2*i,2*i)=dexp(cmplx(0.0d0,-theta(i)))
!end do

do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    tau(i,j)=exp(cmplx(0.0d0,theta((i+1)/2)))
   else
    tau(i,j)=exp(cmplx(0.0d0,-theta(i/2)))
   end if
  else 
   tau(i,j)=0.0d0
  end if
 end do
end do


!call WRCRN('tau',tau,8,8,8)
SL=matmul(tau,S)
!call WRCRN('SL',SL,8,8,8)
 !次にexp-iphiを行列積する
SLd=SL

phi(1)=0.0d0
phi(2)=phid(1)
phi(3)=phid(2)
phi(4)=phid(3)
do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    expmiphi(i,j)=exp(cmplx(0.0d0,-phi((i+1)/2)))
   else
    expmiphi(i,j)=exp(cmplx(0.0d0,-phi(i/2)))
   end if
  else 
   expmiphi(i,j)=0.0d0
  end if
 end do
end do

!!('expmiphi',test,8,8,8)
SL=matmul(expmiphi,SLd)
 !gを施す
 SLd=SL
do i=1,8
 do j=1,4
  SL(2*j-1,i)=-SLd(2*j,i)
  SL(2*j,i)=SLd(2*j-1,i)
 end do
end do




!tau*を施す
SLd=SL
SL=matmul(dconjg(tau),SLd)
!S*を施す
SLd=SL
SL=matmul(dconjg(S),SLd)
!tau*をかける
SLd=SL
SL=matmul(dconjg(tau),SLd)
!expiphiを行列積する
SLd=SL
do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    expiphi(i,j)=exp(cmplx(0.0d0,phi((i+1)/2)))
   else
    expiphi(i,j)=exp(cmplx(0.0d0,phi(i/2)))
   end if
  else 
   expiphi(i,j)=0.0d0
  end if
 end do
end do

!!('expiphi',test,8,8,8)
SL=matmul(expiphi,SLd)
!call WRCRN('SL',SL,8,8,8)
 !SLにgダガーを施す
 SLd=SL
 do i=1,8
  do j=1,4
 SL(2*j-1,i)=SLd(2*j,i)
  SL(2*j,i)=-SLd(2*j-1,i)
 end do
end do


 
!最後にtauをかける
SLd=SL
SL=matmul(tau,SLd)
!call WRCRN('SL',SL,8,8,8)

!!('test3',test,8,8,8)
!ビーネッカー方程式を解く(固有値だけ)
call EVLCG(SL,EVAL)
!call WRCRN('EVAL',EVAL,8,1,8)
do i=1,8
kikaku(i)=dsqrt(dreal(EVAL(i))**2+dimag(EVAL(i))**2)
EVALa(i)=EVAL(i)
if (kikaku(i)<0.001d0) then
EVAL(i)=0.001d0*EVALa(i)
else
EVAL(i)=EVALa(i)/kikaku(i)
end if
end do

do i=1,8
sig(i)=dsign(1.0d0,DIMAG(EVAL(i)))
EVALD(i)=sig(i)*dsqrt((dreal(EVAL(i))+1)/2)
end do
!!!!('EVALD',EVALD,8,1,8)
EVALDRE=dreal(EVALD)
do i=1,7
 do j=i+1,8
 if(EVALDRE(i)>=EVALDRE(j)) then
 else
 y=EVALDRE(i)
 EVALDRE(i)=EVALDRE(j)
 EVALDRE(j)=y
 end if
 end do
end do
wild=dabs(EVALDRE(4))
end function andreev
function andreev2(phid,S,thetax) result(roa)
implicit none
complex(8),intent(in) :: S(8,8)
real(8),intent(in) :: phid(3),thetax
real(8) roa,kikaku(8)
integer  i,j,n,p,q,m,nm,mh
complex(8),dimension(4) :: cnumber,v1
real(8) A_d(4), B_d(4),IP(4),sig(8),sha(2)
real(8) P_so,y,z,phi1,phi2,phi3,EVALDRE(8),ER(4),EI(4)
real(8),parameter :: pi =3.14159265359
complex(8) check
complex(8),dimension(4) :: l,a,b,ad,bd
complex(8),dimension(8,8) :: c
real(8),dimension(4) ::ar,ai,br,bi,clr,cli,phi,theta
real(8),dimension(4,4) :: d_,e_
real(8),dimension(8,8) :: x,xt !tはtildeの略
complex(8),dimension(8,8) :: ct, cd !dはdashの略
complex(8),dimension(8,8) :: cdt
complex(8),dimension(8,8) :: U,Udag !dagはdaggerの略
complex(8),dimension(8,8) :: lambda,test,Sd
complex(8),dimension(8) :: EVAL,EVALD,EDA,EVALa
complex(8),dimension(8,8) :: EVEC,SL,SLd,expmiphi,expiphi,tau

!SLを構成
!tauを構成し、かける

do i=1,4
theta(i)=thetax
end do


do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    tau(i,j)=exp(cmplx(0.0d0,theta((i+1)/2)))
   else
    tau(i,j)=exp(cmplx(0.0d0,-theta(i/2)))
   end if
  else 
   tau(i,j)=0.0d0
  end if
 end do
end do


SL=matmul(tau,S)
 !次にexp-iphiを行列積する
SLd=SL

phi(1)=0.0d0
phi(2)=phid(1)
phi(3)=phid(2)
phi(4)=phid(3)
do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    expmiphi(i,j)=exp(cmplx(0.0d0,-phi((i+1)/2)))
   else
    expmiphi(i,j)=exp(cmplx(0.0d0,-phi(i/2)))
   end if
  else 
   expmiphi(i,j)=0.0d0
  end if
 end do
end do

!('expmiphi',test,8,8,8)
SL=matmul(expmiphi,SLd)
 !gを施す
 SLd=SL
do i=1,8
 do j=1,4
  SL(2*j-1,i)=-SLd(2*j,i)
  SL(2*j,i)=SLd(2*j-1,i)
 end do
end do




!tau*を施す
SLd=SL
SL=matmul(dconjg(tau),SLd)
!S*を施す
SLd=SL
SL=matmul(dconjg(S),SLd)
!tau*をかける
SLd=SL
SL=matmul(dconjg(tau),SLd)
!expiphiを行列積する
SLd=SL
do i=1,8
 do j=1,8
  if (j==i) then
   if (mod(j,2)==1) then
    expiphi(i,j)=exp(cmplx(0.0d0,phi((i+1)/2)))
   else
    expiphi(i,j)=exp(cmplx(0.0d0,phi(i/2)))
   end if
  else 
   expiphi(i,j)=0.0d0
  end if
 end do
end do

!('expiphi',test,8,8,8)
SL=matmul(expiphi,SLd)

 !SLにgダガーを施す
 SLd=SL
 do i=1,8
  do j=1,4
 SL(2*j-1,i)=SLd(2*j,i)
  SL(2*j,i)=-SLd(2*j-1,i)
 end do
end do


 
!最後にtauをかける
SLd=SL
SL=matmul(tau,SLd)
!('SL',SL,8,8,8)

!('test3',test,8,8,8)
!ビーネッカー方程式を解く(固有値だけ)
call EVLCG(SL,EVAL)
!!!('EVAL',EVAL,8,1,8)
do i=1,8
kikaku(i)=dsqrt(dreal(EVAL(i))**2+dimag(EVAL(i))**2)
EVALa(i)=EVAL(i)
if (kikaku(i)<0.001d0) then
EVAL(i)=0.001d0*EVALa(i)
else
EVAL(i)=EVALa(i)/kikaku(i)
end if
end do

do i=1,8
sig(i)=dsign(1.0d0,DIMAG(EVAL(i)))
EVALD(i)=sig(i)*dsqrt((dreal(EVAL(i))+1)/2)
end do
!!!!('EVALD',EVALD,8,1,8)
EVALDRE=dreal(EVALD)
do i=1,7
 do j=i+1,8
 if(EVALDRE(i)<EVALDRE(j)) then
 y=EVALDRE(i)
 EVALDRE(i)=EVALDRE(j)
 EVALDRE(j)=y
 end if
 end do
end do
roa=dabs(EVALDRE(3)-EVALDRE(4))
end function andreev2


function gs(a,n)  result(e)
 implicit none
 integer,intent(in):: n
complex(8),intent(in)::a(n,n)
complex(8) e(n,n),dotp
 integer k,j,i
  e(1:n,1)=normal_vec2(a(1:n,1:1),n)
  do k=2,n
  if (mod(k,2)==1) then
   e(1:n,k)=a(1:n,k)
   do j=1,k-1
    dotp=dot_product1(dconjg(e(1:n,j)),a(1:n,k))
    e(1:n,k)=e(1:n,k)-dotp*e(1:n,j)
      
   end do
   e(1:n,k)=normal_vec2(e(1:n,k:k),n)
   else
    do i=1,n
     if(mod(i,2)==1) then
      e(i,k)=-dconjg(e(i+1,k-1))
     else
      e(i,k)=dconjg(e(i-1,k-1))
     endif
    end do
   end if
  end do
end function gs

function normal_vec2(v,n) result(nv)
 implicit none
 integer,intent(in) :: n
complex(8),intent(in) :: v(n)
complex(8) nv(n),v1
 v1=sqrt(dot_product1(dconjg(v),v))
 if(abs(v1)<10**(-3)) then
   nv(:)=0.0d0
 else
   nv(:)=v(:)/v1
 endif
end function normal_vec2

function gsr(a,n)  result(e)
 implicit none
 integer,intent(in):: n
 real(8),intent(in)::a(n,n)
 real(8) e(n,n),dotp
 integer k,j
  e(1:n,1)=normal_vec2r(a(1:n,1:1),n)
  do k=2,n
   if (mod(k,2)==1) then
    e(1:n,k)=a(1:n,k)
    do j=1,k-1
     dotp=dot_product2(e(1:n,j),a(1:n,k))
     e(1:n,k)=e(1:n,k)-dotp*e(1:n,j)
   
    end do
    e(1:n,k)=normal_vec2r(e(1:n,k:k),n)
   else
    do i=1,n
     if(mod(i,2)==1) then
      e(i,k)=-(e(i+1,k-1))
     else
      e(i,k)=(e(i-1,k-1))
     endif
    end do
   end if
  end do
end function gsr

function normal_vec2r(v,n) result(nv)
 implicit none
 integer,intent(in) :: n
 real(8),intent(in) :: v(n)
 real(8) nv(n),v1
 v1=sqrt(dot_product2(v,v))
 if(abs(v1)<10**(-3)) then
   nv(:)=0.0d0
 else
   nv(:)=v(:)/v1
 endif
end function normal_vec2r

function dot_product1(v,w) result(dotp)
implicit none
integer  n,i
complex(8),intent(in) :: v(:),w(:)
complex(8) dotp
n=size(v(:))
dotp=0.0d0
do i=1,n
dotp=dotp+v(i)*w(i)
end do
end function dot_product1

function dot_product2(v,w) result(dotp)
implicit none
integer  n,i
real(8),intent(in) :: v(:),w(:)
real(8) dotp
n=size(v(:))
dotp=0.0d0
do i=1,n
dotp=dotp+v(i)*w(i)
end do
end function dot_product2

function chirasu1(v) result(w)
implicit none
integer i
real(8),intent(in) :: v(:)
real(8) w(size(v))
do i=1,size(v)
    w(i)=2*v(i)-1
end do  
end function chirasu1

function chirasu2(v) result(w)
implicit none
integer i,j
real(8),intent(in) :: v(:,:)
real(8) w(size(v,1),size(v,2)),sha(2)
sha=shape(v)
do i=1,sha(1)
  do j=1,sha(2)
    w(i,j)=2*v(i,j)-1
  end do
end do  
end function chirasu2

function chirasu3(v) result(w)
implicit none
integer i,j,k
real(8),intent(in) :: v(:,:,:)
real(8) w(size(v,1),size(v,2),size(v,3)),sha(3)
sha=shape(v)
do i=1,sha(1)
  do j=1,sha(2)
   do k=1,sha(3)
    w(i,j,k)=2*v(i,j,k)-1
   end do
  end do
end do  
end function chirasu3
end program wile
