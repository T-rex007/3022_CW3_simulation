{ MAGNET_COIL.PDE  
  

 
  According to Maxwell's equations, 
        curl H = J 
        div B = 0 
        B = mu*H or H = B/mu
 
  where B is the manetic flux density
        H is the magnetic field strength
        J is the electric current density
  and  mu is the magnetic permeability of the material.
 
  The magnetic vector potential A is related to B by 
        B = curl A 
  therefore 
        curl( (1/mu)*curl A ) = J
    if rmu = 1/mu
    then 
          curl( (rmu)*curl A ) = J    
 
  This equation is usually supplmented with the Coulomb Gauge condition 
        div A = 0.
 
  In the axisymmetric case, the current is assumed to flow only in the
  azimuthal direction, and only the azimuthal component of the vector
  potential is present.  Henceforth, we will simply refer to this component as A.
 
  The Coulomb Gauge is identically satisfied, and the PDE to be solved in this
  model takes the form 
        curl((1/mu)*curl (A)) = J(x,y)     in the domain
                           A  = g(x,y)     on the boundary.
 
  The magnetic induction B takes the simple form 
        B = (-dz(A), 0, dr(A)+A/r)
 
  and the magnetic field is given by 
        H = (-dz(A)/mu, 0, (dr(A)+A/r)/mu)
 
  Expanding the equation in cylindrical geometry results in the final equation, 
        dz(dz(A)/mu) + dr((dr(A)+A/r)/mu) = -J
 
  The interpretation of the natural boundary condition becomes 
        Natural(A) = n X H
 
  where n is the outward surface-normal unit vector.
 
  Across boundaries between regions of different material properties, the
  continuity of (n X H) assumed by the Galerkin solver implies that the
  tangential component of H is continuous, as required by the physics.
 
 
  In this simple test problem, we consider a circular coil whose axis of
  rotation lies along the X-axis. We bound the coil by a distant spherical
  surface at which we specify a boundary condition (n X H) = 0.
  At the axis, we use a Dirichlet boundary condition A=0.
}
 
title 'AXI-SYMMETRIC MAGNETIC FIELD'
 
coordinates
    { Cylindrical coordinates, with cylinder axis along Cartesian X direction }
    xcylinder(Z,R)  
 
variables
    Aphi        { the azimuthal component of the  vector potential }
 
definitions
   Lx =0.1
   Ly = 0.1
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
    muu_skin = pi *4e-7 
    muu_muscle=pi *4e-7 
    muu_fat = pi *4e-7 
    muu_cb =  pi *4e-7 
    muu_canb =  pi *4e-7 
    muu_in = pi *4e-7 
    muu_out = pi *4e-7 
    muu_di = pi *4e-7 
    scale = 1
    
     wire_rad = (8.44e-3) * scale
     cop_out = (6.3e-3 - 4.7e-3)* scale                                                               !thickness of outside component
     di_s_d = (4.7e-3 - 0.7e-3) * scale                                                                 !thickness of dielectric
     cop_in = (1.8e-3) * scale
     hole = wire_rad - cop_in - di_s_d -cop_out
  
    current = 2.14
    Bz = dr(r*Aphi)/r                                                                                       ! The Z-component of the Flux density
    Br = dz(Aphi)                                                                                            ! The R-component of the Flux density Note: the phi-component
    Bmag = sqrt(Bz^2 + Br^2)                                                                         ! The Magnitude of the flux density
    H = curl(Aphi)/mu                                                                                      ! The magnetic field strength
    Hmag = magnitude(H)                                                                               ! The magnitude magnetic field strength
    B = curl(Aphi)                                                                                            ! The magnetic flux density
    !Bmag = magnitude(B)
    !J1 = curl(H)
initial values
    Aphi = 2            { unimportant unless mu varies with H }
 
equations
    { FlexPDE expands CURL in proper coordinates }
    Aphi : curl(rmu*curl(Aphi)) = J 
 
boundaries
    region 1 'Domain'
     mu = mu0
       START 'ring' ( -Lx , 0)
      value(Aphi) = 0       { specify A=0 along axis }
        line to (Lx,0)
      natural(Aphi) = 0     { H<dot>n = 0 on distant sphere }
        arc(center=0,0) angle 180 to close         !forms semicircle on Z axis with radius of 0.1 to fully encapsulate the ring
    
    region 2 'Skin'
    	mu = muu_skin
        start (0, ar)
        line to (sth,  ar) to (sth,  0) to (0, 0) to close
        
   region 3 'Fat'
    	mu = muu_fat
        start (sth, ar)
        line to (sth + fth,  ar) to (sth + fth,  0) to (sth, 0) to close
        
    region 4 'Muscle'
    	mu = muu_muscle 
        start (sth + fth, ar)
        line to (sth + fth + mth,  ar) to (sth + fth + mth,  0) to (sth + fth, 0) to close
        
    region 5 'Cortical Bone'
    	mu = muu_cb
        start (sth + fth + mth, ar)
        line to (sth + fth + mth + dcb,  ar) to (sth + fth + mth + dcb,  0) to (sth + fth + mth, 0) to close
        
    region 6 'Cancellous Bone'
    	mu = muu_canb
        start (sth + fth + mth +dcb, ar)
        line to (sth + fth + mth + dcb + dcanb,  ar) to (sth + fth + mth + dcb + dcanb,  0) to (sth + fth + mth + dcb, 0) to close
        
    {region 7 'Hole'
        mu = mu0
       START 'ring' (pen_len , hole)
           line to (pen_len, 0)
           line to (-wire_length , 0)
           line to (-wire_length ,  hole)
           line to close} 
           
    region 8 'Terminal'
        mu = muu_in
        J = current/(pi*((hole + cop_in)^2 - hole^2))
       START 'ring' (pen_len , hole + cop_in)
           line to (pen_len, hole)
           line to (-wire_length , hole)
           line to (-wire_length ,  hole + cop_in)
           line to close 
        
    region 9 'Dielectric'
        mu = muu_di
       START 'ring' (pen_len , wire_rad )
           line to (pen_len, cop_in)
           line to (-wire_length , cop_in)
           line to (-wire_length ,  wire_rad)
           line to close 
           
    region 10 'Ground'
        mu = muu_out
        J = -current/(pi*((hole + cop_in + cop_out)^2 - (hole + cop_in)^2))
       START 'ring' (pen_len, wire_rad)
           line to (pen_len, wire_rad - cop_out)
           line to (0 , wire_rad - cop_out)
           line to (0 ,  wire_rad)
           line to close 
           
monitors
    contour(Bz) zoom(-2,0,4,4) as 'FLUX DENSITY B'
    contour(Aphi) as 'Potential'
 
plots
    grid(z,r)
    contour(Bz)  as 'FLUX DENSITY B'
    contour(Bz) zoom(-2,0,4,4)  as 'FLUX DENSITY B'
    contour(Bmag) as 'Magnitude of the Flux Density (T)' PNG
    contour(Hmag) as  'Magnitude of Field intensity (A/m)' PNG
    elevation(Aphi,dr(Aphi),Aphi/r,dr(Aphi)+Aphi/r,Aphi+r*dr(Aphi)) 
        from (0,0) to (0,1) as 'Bz'
    vector(dr(Aphi)+Aphi/r,-dz(Aphi)) as 'FLUX DENSITY B'
    vector(dr(Aphi)+Aphi/r,-dz(Aphi)) zoom(-2,0,4,4) as 'FLUX DENSITY B'
    vector(H) as 'Vector Plot of magnetic field intensity H (A/m)' PNG
    vector(B) as 'Vector plot of the Manetic Flux Density B (T)' PNG
    vector(Aphi)  as 'VECTOR PLOT OF MAGNETIC POTENTIAL' PNG
    vector(J) zoom(0, 0, pen_len, wire_rad)
    contour(Aphi)  as 'MAGNETIC POTENTIAL' PNG
    contour(Aphi) zoom(-2,0,4,4)  as 'MAGNETIC POTENTIAL'
    
    surface(Aphi)  as 'MAGNETIC POTENTIAL'  viewpoint (-1,1,30)
    
{Summary
    !report (GlobalMax(Bmag, 1)) as 'Maximum flux density  (Wb/m^2)'
    !report (GlobalMax(Bmag, 2)) as 'Maximum flux density in the wire (Wb/m^2)'
    report (GlobalMax(Bmag, 3)) as 'Maximum flux density in the skin (Wb/m^2)'
    report (GlobalMax(Bmag, 4)) as 'Maximum flux density in the fat (Wb/m^2)'
    report (GlobalMax(Bmag, 5)) as 'Maximum flux density in the muscle (Wb/m^2)'
    
    !report (GlobalMax(Hmag, 1)) as 'Maximum Magnetic Field intensity (A/m)'
    !report (GlobalMax(Hmag, 2)) as 'Maximum Magnetic Field intensity in the wire (A/m)'
    report (GlobalMax(Hmag, 3)) as 'Maximum Magnetic Field intensity in the skin (A/m)'
    report (GlobalMax(Hmag, 4)) as 'Maximum Magnetic Field intensity in the fat (A/m)'
    report (GlobalMax(Hmag, 5)) as 'Maximum Magnetic Field intensity in the muscle (A/m)'
}end "CIEPcbP1p8UrD0Ykk9rG4yeq+J6b384kTZRnemmylDzh4NSeLTdnTMCRkOBMoS23tpclJAzvcbEpb3QAb6K7Qovjmwqs32XquSq7Zn/zuZvD0C8taSgsHPSQzqwzfwJ8wpGqtCfyHAV+cyJ02PPcvZ8Sc5HVbAp95qxUD3S3CqH"
