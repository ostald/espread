using LinearAlgebra
include("constants.jl")


function ode_boris_mover_mfp(n_mfp, r0v0, q, m, Bin, sigma, densityf; wMax = [], OPS = [])
  # Boris Mover adapted
  #   - independent variable is number of mean free paths of a particle instead of time
  #   - include ionization cross sections of N2, O2, O and densities for mean free path calculations
  #   - removed electric field

  dOPS = (nPerGyro = 20,
          #event_fcn=[], 
          #StopAfterEvent = 1,
          #i_save = 191 #save every 191st iteration 
          );

  """
  if OPS != [] #if nargin > 7
    dOPS = (; dOPS..., OPS...) #merge_structs(dOPS,OPS);
  end
  """

  nPerGyro = dOPS.nPerGyro;
  #i_save = dOPS.i_save

  """
  if isempty(dOPS.event_fcn)
    hasEventFcn = false;
    rvt_e = [];
  else
    hasEventFcn = true;
    event_fcn = dOPS.event_fcn;
    iEvent = 1;
    r_prev = r0v0[1:3];
    v_prev = r0v0[4:6];
    StopAfterEvent = dOPS.StopAfterEvent;
    rvt_e = [];
  end
  """

  if isempty(wMax) #if nargin < 7 || isempty(wMax)
    wMax = 0;
  end

  ## Start assigning output variables:
  #rvt = zeros((7, 1000)) #rv[6,numel(t)] = 0;
  #rvt[:,1] = [r0v0; 0]
  
  # Initialize the particle position and velocity:
  r = r0v0[1:3];
  v = r0v0[4:6];
  #r = rvt[1:3,1];
  #v = rvt[4:6,1];
  Ekin = E_ev(norm(v))
  cross_sections = sigma(Ekin)


  # Starting time:
  #total_path = 0
  t_running = 0;
  #i_next = 2;
  td_factor = 1 ./ nPerGyro;
  #t = 0:1:600
  mfp_running = 0

  while mfp_running < n_mfp
    B = Bin(r);

    ## Update the increment in time
    #  to ensure nice trajectories along the track use 1 20th (or
    #  what value we give the parameter: nPerGyro) of a gyro-period,
    #  or the time to the next requested time for output.
    #  
    #  Start with the gyro-frequency:
    w_c = abs(q*norm(B)/m);
    Dt = td_factor/(max(w_c,wMax)/2/pi);#t[i_next]-t_running);
    
    ## Then the Boris-scheme:
    q_prime = Dt*q/(2*m);
    h = q_prime*B;
    s = 2*h/(1+norm(h)^2);
    u = v;
    u_prime = u + cross( u + cross( u, h ), s);
    
    v_halfstep = u_prime;
    r = r + v_halfstep*Dt;
    v = v_halfstep;
    
    """
    if hasEventFcn
      EventHappened = event_fcn(t,r,v,E,B,r_prev,v_prev);
      if EventHappened
        rvt_e[iEvent].B = B;
        rvt_e[iEvent].E = E;
        rvt_e[iEvent].r = r;
        rvt_e[iEvent].v = v;
        rvt_e[iEvent].r_prev = r_prev;
        rvt_e[iEvent].v_prev = v_prev;
        rvt_e[iEvent].t = t_running;
        iEvent = iEvent + 1;
        if StopAfterEvent
          rvt = rvt[:,1:[i_next-1]];
          return [rvt, rvt_e]
        end
      end
      r_prev = r;
      v_prev = v;
    end
    """
    t_running = t_running + Dt;

    #assume path segments are short enough to be approximated by
    #total_path = total_path + norm(v)*Dt;

    #local mean freep path:
    altitude = norm(r) - c.re
    if isnan(densityf(altitude)[1])
      println("Altitude [m]: " * string(altitude))
      error("Density is nan. Check altitude range")
    end
    lam = sum(cross_sections .* densityf(altitude))
    
    #mfp = 1/sum(sigma(Ekin) .* ns(altitude))


    #path travelled in fractions of local mfp:
    dr_mfp = norm(v)*Dt * lam
    #f_mfp = norm(v)*Dt / mfp

    #running sum of mfp
    mfp_running = mfp_running + dr_mfp
    #println(string(i_next))
    #println(string(t_running))
    #println(string(mfp_running))
    # create t_save parameter
    
    """
    if i_next % i_save == 0
      ind = Int(i_next/i_save)
      rvt[:,ind] = [r[:];v[:];t_running]
      println(string(cross_sections))
      println(string(dr_mfp))
      println(string(Dt))
      println(string(lam))
      println(string(v))
      println(string(i_next))
      println(string(t_running))
      println(string(mfp_running))
      println(" ")
      if size(rvt, 2) == ind
        rvt = hcat(rvt, zeros((7, 1000)))
      end
    end
    """
    #i_next = i_next+1;
  end
  return r, v, t_running
  #last_ind = ceil(Int, i_next/i_save)
  #rvt[:,last_ind] = [r[:];v[:];t_running]
  #println(string(last_ind))
  #return [rvt[:, 1:last_ind], rvt_e, t_running]
end




function ode_boris_mover(t, r0v0, q, m, Ein, Bin; wMax = [], OPS = [])
    # ODE_BORIS_MOVER - charged particle equation-of-motion-solver in B and E-fields
    #   The ode_boris_mover integrates the equations of motion in
    #   magnetic and electrical fields with the Boris-mover scheme
    #   [e.g. 1].
    # 
    # Calling:
    #   rv, re = ode_boris_mover(t,r0v0,q,m,Ein,Bin[,wMax][,OPS])
    # Input:
    #   t    - time (s), double array [1 x n_t] for desired output
    #   r0v0 - initial particle position, double array [6 x 1] with the
    #          first 3 components for the position (m) and the last 3
    #          positions for the velocity (m/s). The velocities will be
    #          treated as the velocity at half the time-step before t0,
    #          that is no special treatment of the velocities at the
    #          initiation.
    #   q    - particle charge (C), scalar double
    #   m    - particle mass (kg), scalar double
    #   Ein  - electrical field (V/m), either a [3 x 1] array for a
    #          static uniform E-field or a function handle to a
    #          function f_E(t,r) returning a [3 x 1] electrical in
    #          position r at time t.
    #   Bin  - magnetic field (T), either a [3 x 1] array for a
    #          static uniform B-field or a function handle to a
    #          function f_B(t,r) returning a [3 x 1] electrical in
    #          position r at time t.
    #   wMax - maximum expected gyro-frequency, used to limit the
    #          time-step of the solver. Optional argument, defaults to
    #          zero, for which the time-step is controlled by the B-field at
    #          current position and the requested number of steps per
    #          gyro-orbit.
    #   OPS  - struct with options controlling the behaviour. Current fields
    #          are: 
    #          .nPerGyro: number of time-steps per gyro-period, defaults to 20
    #          .event_fcn: function-handle to events-function (change of
    #              direction in magnetic mirror, crossing some boundary, etc).
    #              The function should return 1 if the event has occurred
    #              between the one time-step and the next and have the call:
    #              event_fcn(t,r,v,E,B,r_prev,v_prev). The event-function will
    #              be called at each time-step of the boris-mover, that is at
    #              the the points the boris-mover steps, not the points that
    #              are requested for the resulting trajectory. This limits the
    #              precision of the event to T_gyro/nPerGyro. Further
    #              refinement by interpolation or focused integration will have
    #              to be done by the user, possibly in repeated calls to
    #              ode_boris_mover with the outputs in rv_e as initial
    #              conditions and a finer time-grid. Default: empty - no events
    #          .StopAfterEvent: true/false flag controlling whether to stop the
    #              integration after an event or continue. Defaults to 1.
    # Output:
    #  rv    - particle phase-space coordinates, double array [6 x n_t].
    #          The first 3 rows are the particle position (m), rows 4
    #          to 6 are the particle velocity (m/s),
    #  rv_e  - struct array with the parameters at times just after and before
    #          events. Fields are: r (position (m), [3 x 1] double array), v
    #          (velocity (m/s), [3 x 1] double array), E (E-field V/m, double
    #          array [3 x 1]), B  (B-field T, double array [3 x 1], these
    #          parameters are for the time-step just after the event, r_prev
    #          and v_prev - position and velocities from the previous
    #          time-step. If no events occurs an empty array is returned.
    #  OPS   - default options-struct, returned when function called without
    #          input arguments.
    # 
    # The equations of the Boris scheme are:
    #   x_{k+1}   = x_k + Dt*v_{k+1/2}
    #   v_{k+1/2} = u' + q' E_k
    #
    # with
    #
    #   u' = u + ( u + ( u x h ) ) x s
    #   u  = v_{k-1/2} + q'*E_k
    #   h  = q'*B_k
    #   s  = 2*h/( 1 + h^2 )
    #
    #   q' = Dt*(q/2m)
    # 
    # Example:
    #   # Physical constants and parameters
    # m_e = 9.1094e-31;        # electron mass (kg)
    # q_e = 1.6022e-19;        # electron charge (C)
    # B = 5e-5;                # Magnetic field strength (T)
    # w_e = (q_e*B/m_e);       # electron gyro-frequency
    # T_gyro = 1/(w_e/(2*pi)); # gyro-period
    # # Initial conditions and time-span
    # v_0    = (2*q_e*1/m_e).^(1/2); # velocity of 1 eV electron (m/s)
    # r_gyro = v_0/w_e;               # electron gyro-radius
    # t_out = LinRange(0,100*T_gyro,10001); # requested time-steps 
    # 
    # # ODE-integration:
    # rv, re = ode_boris_mover(t_out,[0;-r_gyro;0;v_0[1];0;0],-q_e,m_e,[0;0;0],[0;0;B]);
    # 
    # # Digestion of solution
    # using Statistics
    # r_c = mean(rv[1:3,:], dims=2); #center of 
    # using GLMakie
    # 
    # ##
    # 
    # f = Figure()
    # ax1 = Axis(f[1, 1],
    #    ylabel = "y-position (m)",
    #    xlabel = "x-position (m)",
    #    title = "X-Y trajectory 100 gyrations",
    #    ) #subplot(2,2,1)
    # lines!(ax1, rv[1,:], rv[2,:])
    # 
    # ax2 = Axis(f[1, 2],
    #    ylabel = "y-velocity (m/s)",
    #    xlabel = "x-velocity (m/s)",
    #    title = "v_x-v_y 100 gyrations",
    #    ) #subplot(2,2,3)
    # lines!(ax2, rv[4,:], rv[5,:])
    # 
    # r_max = maximum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
    # r_min = minimum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
    # r_g =   mean(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
    # ax3 = Axis(f[2, 1],
    #    ylabel = "(r_max - r_min)/r_g",
    #    xlabel = "time (s)",
    #    title = "relative variation of gyro-radius: $(round((r_max-r_min)/r_gyro, digits=6))",
    #    ) #subplot(2,2,3)
    # lines!(ax3, t_out,((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]) .^2 ) .^.5); #???
    # 
    # K_max = maximum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e);
    # K_min = minimum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e);
    # ax4 = Axis(f[2, 2],
    #    ylabel = "(K_{max} - K_{min})/K_0",
    #    xlabel = "time (s)",
    #    title = "relative variation of Kinetic energy: $(round((K_max-K_min)/1, digits = 16))",
    #    ) #subplot(2,2,4)
    # lines!(ax4, t_out,(rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e)
    # 
    #   The current implementation adjusts the time-increment to 1/20
    #   of the gyro-period at the current particle position. That might
    #   lead to errors where the B-field varies abruptly over distances
    #   smaller than |1/r_gyro|. The mover conserves kinetic energy
    #   for particles in a uniform B-field to machine precision, and
    #   keeps relative variation of gyro-radius to ~4e-4.
    
    #  Copyright � Bjorn Gustavsson 20190118, bjorn.gustavsson@uit.no
    #  This is free software, licensed under GNU GPL version 2 or later

    #defining default parameters:

    dOPS = (nPerGyro = 20,
            event_fcn=[], 
            StopAfterEvent = 1,
            );

    """
    if OPS != [] #if nargin > 7
      dOPS = (; dOPS..., OPS...) #merge_structs(dOPS,OPS);
    end
    """

    nPerGyro = dOPS.nPerGyro;
    
    if isempty(dOPS.event_fcn)
      hasEventFcn = false;
      rv_e = [];
    else
      hasEventFcn = true;
      event_fcn = dOPS.event_fcn;
      iEvent = 1;
      r_prev = r0v0[1:3];
      v_prev = r0v0[4:6];
      StopAfterEvent = dOPS.StopAfterEvent;
      rv_e = [];
    end

    if isempty(wMax) #if nargin < 7 || isempty(wMax)
      wMax = 0;
    end

    ## Start assigning output variables:
    rv = zeros((6, length(t))) #rv[6,numel(t)] = 0;
    rv[:,1] = r0v0;
    
    # Initialize the particle position and velocity:
    r = rv[1:3,1];
    v = rv[4:6,1];
    # Starting time:
    t_running = t[1];
    i_next = 2;
    td_factor = 1 ./ nPerGyro;
    while t_running < t[end]

    # Electric field is 0 for my application, so 
      
      # Calculate the E and B-fields at current positions:
      if isempty(Ein) # "Special case for when the E-field depends on
                      # the local B-field (case: plasma-spheric
                      # corotation) - then the Bin function has to
                      # return both B and E fields
        B, E = Bin(t_running,r);
      else
        if isa(Ein, Function) #if isa(Ein,'function_handle')
          # Then we have a function for the E-field that might depend
          # on both time and position
          E = Ein(t_running,r);
        else
          # we have a constant E-field specified as a 3-by-1 array
          E = Ein;
        end
        if isa(Bin, Function) #if isa(Bin,'function_handle')
          # Then we have a function for the magnetic field that might
          # depend on both time and position.
          B = Bin(t_running,r);
        else
          B = Bin;
        end
      end

      ## Update the increment in time
      #  to ensure nice trajectories along the track use 1 20th (or
      #  what value we give the parameter: nPerGyro) of a gyro-period,
      #  or the time to the next requested time for output.
      #  
      #  Start with the gyro-frequency:
      w_c = abs(q*norm(B)/m);
      Dt = min(td_factor/(max(w_c,wMax)/2/pi), t[i_next]-t_running);
      
      ## Then the Boris-scheme:
      q_prime = Dt*q/(2*m);
      h = q_prime*B;
      s = 2*h/(1+norm(h)^2);
      u  = v + q_prime*E;
      u_prime = u + cross( u + cross( u, h ), s);
      
      v_halfstep = u_prime + q_prime*E;
      r = r + v_halfstep*Dt;
      v = v_halfstep;
      if hasEventFcn
        EventHappened = event_fcn(t,r,v,E,B,r_prev,v_prev);
        if EventHappened
          rv_e[iEvent].B = B;
          rv_e[iEvent].E = E;
          rv_e[iEvent].r = r;
          rv_e[iEvent].v = v;
          rv_e[iEvent].r_prev = r_prev;
          rv_e[iEvent].v_prev = v_prev;
          rv_e[iEvent].t = t_running;
          iEvent = iEvent + 1;
          if StopAfterEvent
            rv = rv[:,1:[i_next-1]];
            return [rv, rv_e]
          end
        end
        r_prev = r;
        v_prev = v;
      end
      t_running = t_running + Dt;
      if t_running == t[i_next]
        rv[:,i_next] = [r[:];v[:]];
        i_next = i_next+1;
      end
      
    end
    
    return [rv,rv_e]
end



# to do:
# test!

"""

struct DotStruct
  properties::Dict{Symbol, Any}
end
Dstruct() = DotStruct(Dict{Symbol, Any}())
Base.getproperty(x::DotStruct, property::Symbol) = getfield(x, :properties)[property]
Base.setproperty!(x::DotStruct, property::Symbol, value) = getfield(x, :properties)[property] = value
Base.propertynames(x::DotStruct) = keys(getfield(x, :properties))

##
n_mfp = 1e-2
E_0 = 50 #EV
v_0    = (2*c.qe*E_0/c.me).^(1/2); 
r0v0 = [0, 0, c.re+500e3, v_0, 0, 0]
Bin(p) = [0, 0, 4.5e-5]

rvt, rvt_e, t_end = ode_boris_mover_mfp(n_mfp, r0v0, -c.qe, c.me, Bin, sigma_tot, densityf)

w_e = (c.qe*norm(Bin([0, 0, 0]))/c.me);       # electron gyro-frequency
r_gyro = v_0/w_e;               # electron gyro-radius

t_out = rvt[7, :]
rv = rvt[1:6, :]
r_c = mean(rv[1:3,:], dims=2)
f = Figure()
ax1 = Axis(f[1, 1],
   ylabel = "y-position (m)",
   xlabel = "x-position (m)",
   title = "X-Y trajectory 100 gyrations",
   ) #subplot(2,2,1)
scatter!(ax1, rv[1,:], rv[2,:])

ax2 = Axis(f[1, 2],
   ylabel = "y-velocity (m/s)",
   xlabel = "x-velocity (m/s)",
   title = "v_x-v_y 100 gyrations",
   ) #subplot(2,2,3)
scatter!(ax2, rv[4,:], rv[5,:])

r_max = maximum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
r_min = minimum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
r_g =   mean(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
ax3 = Axis(f[2, 1],
   ylabel = "(r_max - r_min)/r_g",
   xlabel = "time (s)",
   title = "relative variation of gyro-radius: _dollar_(round((r_max-r_min)/r_gyro, digits=6))",
   ) #subplot(2,2,3)
lines!(ax3, t_out,((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]) .^2 ) .^.5); #???

K_max = maximum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * c.me / 2 /c.qe);
K_min = minimum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * c.me / 2 /c.qe);
ax4 = Axis(f[2, 2],
   ylabel = "(K_{max} - K_{min})/K_0",
   xlabel = "time (s)",
   title = "relative variation of Kinetic energy: _dollar_round((K_max-K_min)/1, digits = 16))",
   ) #subplot(2,2,4)
lines!(ax4, t_out,(rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * c.me / 2 /c.qe);
display(GLMakie.Screen(), f)

##


m_e = 9.1094e-31;        # electron mass (kg)
q_e = 1.6022e-19;        # electron charge (C)
B = 5e-5;                # Magnetic field strength (T)
w_e = (q_e*B/m_e);       # electron gyro-frequency
T_gyro = 1/(w_e/(2*pi)); # gyro-period
# Initial conditions and time-span
E_0 = 50 #EV
v_0    = (2*q_e*E_0/m_e).^(1/2); # velocity of 1 eV electron (m/s)
r_gyro = v_0/w_e;               # electron gyro-radius
t_out = LinRange(0,100e3*T_gyro,10001); # requested time-steps 

# ODE-integration:
rv, re = ode_boris_mover(t_out,[0;-r_gyro;0;v_0[1];0;0],-q_e,m_e,[0;0;0],[0;0;B]);

# Digestion of solution
using Statistics
r_c = mean(rv[1:3,:], dims=2); #center of 
using GLMakie


f = Figure()
ax1 = Axis(f[1, 1],
   ylabel = "y-position (m)",
   xlabel = "x-position (m)",
   title = "X-Y trajectory 100 gyrations",
   ) #subplot(2,2,1)
lines!(ax1, rv[1,:], rv[2,:])

ax2 = Axis(f[1, 2],
   ylabel = "y-velocity (m/s)",
   xlabel = "x-velocity (m/s)",
   title = "v_x-v_y 100 gyrations",
   ) #subplot(2,2,3)
lines!(ax2, rv[4,:], rv[5,:])

r_max = maximum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
r_min = minimum(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
r_g =   mean(((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]).^2) .^.5);
ax3 = Axis(f[2, 1],
   ylabel = "(r_max - r_min)/r_g",
   xlabel = "time (s)",
   title = "relative variation of gyro-radius: _dollar_(round((r_max-r_min)/r_gyro, digits=6))",
   ) #subplot(2,2,3)
lines!(ax3, t_out,((rv[1,:] .- r_c[1]) .^2 + (rv[2,:] .- r_c[2]) .^2 ) .^.5); #???

K_max = maximum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e);
K_min = minimum((rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e);
ax4 = Axis(f[2, 2],
   ylabel = "(K_{max} - K_{min})/K_0",
   xlabel = "time (s)",
   title = "relative variation of Kinetic energy: _dollar_(round((K_max-K_min)/1, digits = 16))",
   ) #subplot(2,2,4)
lines!(ax4, t_out,(rv[4,:] .^2 + rv[5,:] .^2 + rv[6,:] .^2) * m_e / 2 /q_e)
display(GLMakie.Screen(), f)
"""