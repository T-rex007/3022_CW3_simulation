TITLE 'Heat flow around an Insulating Canister'
COORDINATES
Cartesian3

variables
    Ax        { the azimuthal component of the  vector potential }
    Ay        { the azimuthal component of the  vector potential }
 
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
  
    current = 2.14
    B = curl(Ax,Ay,0)                                                                                         { magnetic flux density }
    Bmag = magnitude(B)                                                                        ! The Magnitude of the flux density
    H = B/mu                                                                                      ! The magnetic field strength
    Hmag = magnitude(H)                                                                               ! The magnitude magnetic field strength

    
initial values
    Ax = 2            { unimportant unless mu varies with H }
    Ay = 2            { unimportant unless mu varies with H }
 
equations
    { FlexPDE expands CURL in proper coordinates }
    Ax: curl(rmu*curl(Ax)) = J
    Ay: curl(rmu*curl(Ay)) = J 
    
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
      LINE TO (Lx,-Ly) TO (Lx,Ly) TO (-Lx,Ly) TO CLOSE
        
LIMITED REGION 2 'tissue' { the embedded blob }
     LAYER 1 
     mu = muu_canb
     eps = eps_canb*eps0
        
     LAYER 2 
     mu = muu_canb
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
     eps = eps_in*eps0
     J = current/(pi*((hole + cop_in)^2 - hole^2))
        
     LAYER 5 
     mu = muu_in
     eps = eps_in*eps0
     J = current/(pi*((hole + cop_in)^2 - hole^2))
        
     LAYER 6 
     mu = muu_in
     eps = eps_in*eps0
     J = current/(pi*((hole + cop_in)^2 - hole^2))
     
     LAYER 7
     mu = muu_in
     eps = eps_in*eps0
     J = current/(pi*((hole + cop_in)^2 - hole^2))
     
     START(cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(hole, 0)
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
     J = -current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
        
     LAYER 5 
     J = -current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
        
     LAYER 6 
     J = -current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
     
     START(cop_out + di_s_d+cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(di_s_d + cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
PLOTS
GRID(y,z) ON x=0
     CONTOUR(Ax^2 + Ay^2) ON z=0
     VECTOR(B) ON z=0
     !vector(Ex, Ey, Ez)
     !ELEVATION(v) FROM (0,-1,0) to (0,1,0) { note 3D coordinates }
END