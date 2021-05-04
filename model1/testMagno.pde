TITLE 'Heat flow around an Insulating Canister'
COORDINATES
Cartesian3

variables
    Ax        { the azimuthal component of the  vector potential }
    Ay
    Az
 
definitions
   !scale = 1
   scale = 1.2
   !scale = 1.4
   !scale = 1.6
   !scale = 1.8
   !scale = 2
   !scale = 2.2
   Lx =0.05*scale
   Ly = 0.05*scale
    mu0 = pi *4e-7            { the permeability }
    mu 
    J0 = 0
    Jx = 0  Jy = 0  Jz = J0
    J = vector(Jx, Jy, Jz)
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
    eps_skin = 8.58e2
    
    muu_muscle=pi *4e-7 
    eps_muscle = 8.26e2
    
    muu_fat = pi *4e-7 
    eps_fat = 4.46e2
    
    muu_cb =  pi *4e-7 
    eps_cb = 1.06e2
    
    muu_canb =  pi *4e-7
    eps_canb = 1.85e2
    
    muu_in = 800*pi *4e-7
    eps_in = 1
    
    muu_out = 800*pi *4e-7 
    eps_out = 1
    
    muu_di = pi *4e-7 
    eps_di = 2.1
    
    {wire_rad = (1e-3) * scale
    cop_out = (0.2e-3)* scale                                                               !thickness of outside component
    di_s_d = (0.4e-3) * scale                                                                 !thickness of dielectric
    cop_in = (0.2e-3) * scale
    hole = wire_rad - cop_in - di_s_d -cop_out}

     wire_rad = (8.44e-3) * scale
     cop_out = (6.3e-3 - 4.7e-3)* scale                                                               !thickness of outside component
     di_s_d = (4.7e-3 - 0.7e-3) * scale                                                                 !thickness of dielectric
     cop_in = (1.8e-3) * scale
     hole = wire_rad - cop_in - di_s_d -cop_out
  
    current = 1E-3
    B = curl(Ax, Ay, Az)                                                                                         { magnetic flux density }
    Bmag = magnitude(B)                                                                        ! The Magnitude of the flux density
    H = B/mu                                                                                      ! The magnetic field strength
    Hmag = magnitude(H)
    
!initial values
   ! Az=2 
   ! Ax=2 
    !Ay=2            { unimportant unless mu varies with H }
 
equations

    Ax: dx( dx(Ax)/mu)+ dy( dy(Ax)/mu)+ dz( dz(Ax)/mu)=-Jx
    Ay: dx( dx(Ay)/mu)+ dy( dy(Ay)/mu)+ dz( dz(Ay)/mu)=-Jy
    Az: dx( dx(Az)/mu)+ dy( dy(Az)/mu)+ dz( dz(Az)/mu)=-Jz
    
EXTRUSION
SURFACE 'Bottom' z=0
      LAYER 'Cancellous Bone'
SURFACE 'Can top' z=dcanb
      LAYER 'Cortical Bone'
SURFACE 'Cor Top' z=dcb + dcanb
      LAYER 'Muscle'
SURFACE 'Cablebot' z=sth + fth + mth + dcb - dcanb-pen_len
      Layer 'Muscle2'
SURFACE 'Mus Top' z=mth + dcb + dcanb
      LAYER 'Fat'
SURFACE 'Fat top' z=fth + mth + dcb + dcanb
      LAYER 'Skin'
SURFACE 'skin Top' z=sth + fth + mth + dcb + dcanb
SURFACE 'Top' z=0.06*scale

BOUNDARIES
REGION 1 'box' 
      eps = eps0
      mu = mu0
      START(-Lx,-Ly)
      value(Ax) = 0  value(Ay) = 0  value(Az) = 0
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
     J0 = -current/(pi*((hole + cop_in)^2 - hole^2))
        
     LAYER 5 
     mu = muu_in
     eps = eps_in*eps0
     J0 = -current/(pi*((hole + cop_in)^2 - hole^2))
        
     LAYER 6 
     mu = muu_in
     eps = eps_in*eps0
     J0 = -current/(pi*((hole + cop_in)^2 - hole^2))
     
     LAYER 7
     mu = muu_in
     eps = eps_in*eps0
     J0 = -current/(pi*((hole + cop_in)^2 - hole^2))
     
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
     J0 = current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
        
     LAYER 5 
     J0 = current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
        
     LAYER 6 
     J0 = current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
     mu = muu_out
     eps = eps_out*eps0
     
     START(cop_out + di_s_d+cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(di_s_d + cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
PLOTS
GRID(y,z) ON x=0
 !CONTOUR(Ax^2 + Ay^2) ON z=0
     VECTOR(J) ON x=0
     vector(B) ON z = sth + fth + mth + dcb - dcanb-pen_len
     !ELEVATION(v) FROM (0,-1,0) to (0,1,0) { note 3D coordinates }
     
Summary
     report globalmax(Bmag, 'Skin')
     report globalmax(Bmag, 'Fat')
     report globalmax(Bmag, 'Muscle2')
     report globalmax(Bmag, 'Muscle')
     report globalmax(Bmag, 'Cortical Bone')
     report globalmax(Bmag, 'Cancellous Bone')
     report ''
     report globalmax(Hmag, 'Skin')
     report globalmax(Hmag, 'Fat')
     report globalmax(Hmag, 'Muscle2')
     report globalmax(Hmag, 'Muscle')
     report globalmax(Hmag, 'Cortical Bone')
     report globalmax(Hmag, 'Cancellous Bone')
END