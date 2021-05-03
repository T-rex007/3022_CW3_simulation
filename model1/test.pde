TITLE 'Heat flow around an Insulating Canister'
COORDINATES
Cartesian3

variables
    v        { the azimuthal component of the  vector potential }
 
definitions
   Lx =0.005
   Ly = 0.005
    mu0 = pi *4e-7            { the permeability }
   ! mu 
    J = 0              { the source defaults to zero }
   ! rmu = 1/ mu
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
    
    muu_in = pi *4e-7
    eps_in = 1
    
    muu_out = pi *4e-7 
    eps_out = 1
    
    muu_di = pi *4e-7 
    eps_di = 2.1
    
    scale = 1


    wire_rad = (1e-3) * scale
    cop_out = (1e-3 - 0.75e-3)* scale                                                               !thickness of outside component
    di_s_d = (0.75e-3 - 0.65e-3) * scale                                                                 !thickness of dielectric
    cop_in = (0.65e-3 -0.5e-3) * scale
    hole = wire_rad - cop_in - di_s_d -cop_out
  
   Ex=-dx(v)  Ey=-dy(v)  Ez=-dz(v)			        ! These are the x and y components of the gradient function: E = - grad(V), in Cartesian coordinates.
   Eabs=sqrt(Ex^2+Ey^2 + Ez^2) 					! You can find the magnitude of the E field vector using Pythagoras' Theorem
   v0 = 1
    
    ! Mass Density 
    skin_md = 1109
    fat_md = 911
    cortical_md = 1908
    cancelous_md = 1178
    muscle_md= 1090
    
   ! Conductivity 
   sig_skin =  3.71e-2
   sig_fat = 4.48e-2
   sig_muscle = 5.48e-1
   sig_cortical =  2.85e-2
   sig_cancel = 9.71e-2
   
   Emag = -grad(v)!Magnitude(Ex,Ey,Ez)
   ! Emag = sqrt(Ex^2 + Ey^2 + Ez^2)
   ! SAR Calculations
   SAR_skin = (0.7071 * sig_skin * Emag^2)/skin_md
   SAR_fat =  (0.7071 * sig_fat * Emag^2)/fat_md
   SAR_muscle = (0.7071 * sig_muscle * Emag^2)/muscle_md
   SAR_cortical = (0.7071 * sig_cortical * Emag^2)/cortical_md
   SAR_cancel =  (0.7071 * sig_cancel * Emag^2)/cancelous_md
   
   !S_skin =  
   !S_fat = 
   !S_muscle =
   !S_cortical = 
   !S_cancel = 
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
      !mu = mu0
      START(-Lx,-Ly)
        natural(v)=0 LINE TO (Lx,-Ly)
        NATURAL(v)=0 LINE TO (Lx,Ly)
        natural(v)=0 LINE TO (-Lx,Ly)
        NATURAL(v)=0 LINE TO CLOSE
        
LIMITED REGION 2 'tissue' { the embedded blob }
     LAYER 1 
     !mu = muu_canb
     eps = eps_canb*eps0
        
     LAYER 2 
     !mu = muu_cb
     eps = eps_cb*eps0
        
     LAYER 3 
     !mu = muu_muscle
     eps = eps_muscle*eps0
    
     LAYER 4 
     !mu = muu_muscle
     eps = eps_muscle*eps0
        
     LAYER 5 
     !mu = muu_fat
     eps = eps_fat*eps0
        
     LAYER 6 
     !mu = muu_skin
     eps = eps_skin*eps0
     
LIMITED REGION 3 'Terminal' { the embedded blob }
     LAYER 4 
     !mu = muu_in
     eps = 1e-15
        
     LAYER 5 
     !mu = muu_in
     eps = 1e-15
        
     LAYER 6 
     !mu = muu_in
     eps = 1e-15
     
     LAYER 7
     !mu = muu_in
     eps = 1e-15
     
     START(cop_in+hole, 0)
     value(v) = v0
     arc(center = 0, 0) angle 360 to close
     
     start(hole, 0)
     value(v) = v0
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 4 'Dielectric_Tissue' { the embedded blob }
     LAYER 4 
     !mu = muu_di
     eps = eps_di*eps0
        
     LAYER 5 
     !mu = muu_di
     eps = eps_di*eps0
        
     LAYER 6 
     !mu = muu_di
     eps = eps_di*eps0
     
     START(di_s_d+cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 5 'Dielectric' { the embedded blob }
     LAYER 7
     !mu = muu_di
     eps = eps_di*eps0
     
     START(wire_rad, 0)
     arc(center = 0, 0) angle 360 to close
     
     start(cop_in+hole, 0)
     arc(center = 0, 0) angle 360 to close
     
LIMITED REGION 6 'Ground' { the embedded blob }
     LAYER 4 
     !mu = muu_out
     eps = 1e-15
        
     LAYER 5 
     !mu = muu_out
     eps = 1e-15
        
     LAYER 6 
     !mu = muu_out
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
     
     SUMMARY
     !report(Emag)
     report(integral(SAR_skin)) as "SAR for the Skin: "
     report(integral(SAR_fat)) as "SAR for the Fat: "
     report(integral(SAR_muscle)) as "SAR for the Muscle: "
     report(integral(SAR_cortical)) as "SAR for the Cortical Bone: "
     report(integral(SAR_cancel)) as "SAR for the Cancelous Bone: "
END