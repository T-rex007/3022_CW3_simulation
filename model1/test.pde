TITLE 'Heat flow around an Insulating Canister'
COORDINATES
Cartesian3

variables
    v        { the azimuthal component of the  vector potential }
 
definitions
   Lx =0.05
   Ly = 0.05
    mu0 = pi *4e-7            { the permeability }
    mu 
    J = 0              { the source defaults to zero }
    rmu = 1/ mu
    nudge = 0.0001
    wire_length = 0.02
    pen_len = 1.5/100
    
    sth =0.0015                                                                                                   ! Thickness of the skin
    ar = 0.05                                                                                                        ! radius of the arm
    fth = 0.0085                                                                                                   ! Thickness of the Fat 
    mth=0.0275                                                                                                   ! Thickness of muscle 
    dcb =0.6/100
    dcanb = ar - sth - mth - fth - dcb 
    ! Radius  of the wire
    eps
    eps0=8.854e-12 							    ! The dielectric permittivity of free space.
    
    muu_skin = pi *4e-7 
    eps_skin = eps0
    muu_muscle=pi *4e-7 
    eps_muscle = eps0
    muu_fat = pi *4e-7 
    eps_fat = eps0
    muu_cb =  pi *4e-7 
    eps_cb = eps0
    muu_canb =  pi *4e-7
    eps_canb = eps0
    muu_in = pi *4e-7
    eps_in = 1
    muu_out = pi *4e-7 
    eps_out = 1
    muu_di = pi *4e-7 
    eps_di = eps0
    scale = 1

     wire_rad = (8.44e-3) * scale
     cop_out = (6.3e-3 - 4.7e-3)* scale                                                               !thickness of outside component
     di_s_d = (4.7e-3 - 0.7e-3) * scale                                                                 !thickness of dielectric
     cop_in = (1.8e-3) * scale
     hole = wire_rad - cop_in - di_s_d -cop_out
  
  Ex=-dx(v)  Ey=-dy(v)  Ez=-dz(v)				    ! These are the x and y components of the gradient function: E = - grad(V), in Cartesian coordinates.
  Eabs=sqrt(Ex^2+Ey^2 + Ez^2) 					! You can find the magnitude of the E field vector using Pythagoras' Theorem
  v0 = 1
  !DEx=eps*Ex    DEy=eps*Ey  DEz = epsEz 		        ! Define the constitutive equation: D = epsilon*E.
  !Dabs=sqrt(DEx^2+DEy^2+ DEz^2)				! Define the magnitude of the D-field vector using Pythagoras' Theorem.
  !Jx = sigma*Ex     Jy = sigma*Ey         ! Define the constitutive equation J = sigma*E
  !Jabs = sqrt(Jx^2+Jy^2)                    ! Define the magnitude of the D-field vector using the Pythagora's Theorem 
 
equations
    { FlexPDE expands CURL in proper coordinates }
     v : div(-eps*grad(v)) = 0		! This is essentially Poisson's equation.
    
EXTRUSION
SURFACE 'Bottom' z=0
      LAYER 'Cancellous Bone'
SURFACE 'Can top' z=dcanb
      LAYER 'Cortical Bone'
SURFACE 'Cor Top' z=dcb + dcanb
      LAYER 'Muscle'
SURFACE 'Cablebot' z=sth + fth + mth + dcb - dcanb-pen_len
SURFACE 'Mus Top' z=mth + dcb + dcanb
      LAYER 'Fat'
SURFACE 'Fat top' z=fth + mth + dcb + dcanb
      LAYER 'Skin'
SURFACE 'skin Top' z=sth + fth + mth + dcb + dcanb
SURFACE 'Top' z=0.08

BOUNDARIES
REGION 1 'box' 
      eps = eps0
      mu = mu0
      START(-Lx,-Ly)
        natural(v)=0 LINE TO (Lx,-Ly)
        NATURAL(v)=0 LINE TO (Lx,Ly)
        natural(v)=0 LINE TO (-Lx,Ly)
        NATURAL(v)=0 LINE TO CLOSE
        
LIMITED REGION 2 'tissue' { the embedded blob }
     LAYER 1 
     mu = muu_canb
     eps = eps_canb*eps0
        
     LAYER 2 
     mu = muu_cb
     eps = eps_cb*eps0
        
     LAYER 3 
     mu = muu_muscle
     eps = eps_muscle*eps0
    
     LAYER 4 
     mu = muu_muscle
     eps = eps_muscle*eps0
        
     LAYER 5 
     mu = muu_fat
     eps = eps_fat*eps0
        
     LAYER 6 
     mu = muu_skin
     eps = eps_skin*eps0
     
LIMITED REGION 3 'Terminal' { the embedded blob }
     LAYER 4 
     mu = muu_in
     eps = 1e-15
        
     LAYER 5 
     mu = muu_in
     eps = 1e-15
        
     LAYER 6 
     mu = muu_in
     eps = 1e-15
     
     LAYER 7
     mu = muu_in
     eps = 1e-15
     
     START(cop_in+hole, 0)
     value(v) = v0
     arc(center = 0, 0) angle 360 to close
     
     start(hole, 0)
     value(v) = v0
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 4 'Dielectric_Tissue' { the embedded blob }
     LAYER 4 
     mu = muu_di
     eps = eps_di*eps0
        
     LAYER 5 
     mu = muu_di
     eps = eps_di*eps0
        
     LAYER 6 
     mu = muu_di
     eps = eps_di*eps0
     
     START(di_s_d+cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 5 'Dielectric' { the embedded blob }
     LAYER 7
     mu = muu_di
     eps = eps_di*eps0
     
     START(wire_rad, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 6 'Ground' { the embedded blob }
     LAYER 4 
     mu = muu_out
     eps = 1e-15
        
     LAYER 5 
     mu = muu_out
     eps = 1e-15
        
     LAYER 6 
     mu = muu_out
     eps = 1e-15
     
     START(cop_out + di_s_d+cop_in+hole, 0)
     value(v) = 0
     arc(center = 0, 0) angle 360 to close
     
     start(di_s_d + cop_in+hole, 0)
     value(v) = 0
     arc(center = 0, 0) angle 360 to close
     
PLOTS
GRID(y,z) ON x=0
     CONTOUR(v) ON x=0
     VECTOR(Ey, Ez) ON x=0
     !vector(Ex, Ey, Ez)
     ELEVATION(v) FROM (0,-1,0) to (0,1,0) { note 3D coordinates }
END