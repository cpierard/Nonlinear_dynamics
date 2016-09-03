#some numerical methods for solving ODEs.

"""
      simple_euler(f, n, x0, dt)
      OUT: array_solutions, array_time

Is the simple Euler's method for solving ODEs. Where `f` is the function of type x' = f(x),
n is the number of steps, `x0` is the initial condition, and `dt` is the time step.

"""
function simple_euler(f::Function, n, x0, dt)

    #Defined two arrays, one for the solution, and another for the times.
    trayectory = [x0]
    time = [0.0]

    #iteration
    for i in 1:n

        x_new = trayectory[i] + dt * f(trayectory[i]) #Rule for iteration.
        push!(trayectory, x_new) #Appends new value to the solutions array.
        push!(time, dt*i) #Appends new time to time array.

    end

    return trayectory, time

end

"""
      improved_euler(f, n, x0, dt)
      OUT: array_solutions, array_time

Similar to simple_euler but with new iteration rule. It is a second order method. `f` is the function of type `x' = f(x)`,
n is the number of steps, `x0` is the initial condition, and `dt` is the time step.

"""
function improved_euler(f::Function, n, x0, dt)

    trayectory = [x0]
    time = [0.0]

    for i in 1:n

        x_aux = trayectory[i] + dt*f(trayectory[i])  #This gives evaluate f in a t+dt step,
        x_new = trayectory[i] + dt*(f(trayectory[i]) + f(x_aux))/2  #then use it two average a better x_new value.
        push!(trayectory, x_new)
        push!(time, dt*i)

    end

    return trayectory, time

end

#Famous Runge-Kutta method. 4 order method

"""
      runge_kutta(f, n, x0, dt)
      OUT: array_solutions, array_time

An ODE integrator with an error of order two. It is a second order method. `f` is the function of type `x' = f(x)`,
n is the number of steps, `x0` is the initial condition, and `dt` is the time step.

"""
function runge_kutta(f::Function, n, x0, dt)

    trayectory = collect(x0)
    time = collect(0.0)

    for i in 1:n

        k1 = f(trayectory[i]) * dt
        k2 = f(trayectory[i] + k1/2) * dt
        k3 = f(trayectory[i] + k2/2) * dt
        k4 = f(trayectory[i] + k3) * dt
        x_new = trayectory[i] + (k1 + 2*k2 + 2*k3 + k4)/6
        push!(trayectory, x_new)
        push!(time, dt*i)

    end

    return trayectory, time
end

"""
      plot_rate(f, xi, xf, interval)
      OUT: plot(xs, fun)

This is a tools that I use for doing a geometric analysis of a diferential equation of the form `x' = f(x)`. `f` is the function, xi is the
initial value of the range, `xf` is the final value, and `interval` is the interval between `xi` and `xf`.
"""

function plot_rate(f::Function, xi, xf, interval) #For a geometric interpretation

    fun = []
    xs = []

    for i in xi:interval:xf

        push!(fun,f(i))
        push!(xs, i)

    end

    plot(xs, fun)
    xlabel("x")
    ylabel("f(x)")
end

"""
      logistic_equation(x0, r)
      OUT: x

The logistic equation as a funtion, where `r` is the rate parameter and `x0` is the initial condition. The function evaluates the `x0`
and throws out the new value for `x`.

"""

function logistic_equation(x0, r)

    x = r*x0*(1-x0)
    return x

end

"""
      logistic_iterator(x0, r, n, yn::Bool)
      OUT:  iteration_serie, evaluations
            println(evaluations[n+1])

            (when yn = true): plot(iteration_serie, evaluations)

This function iterates the logistic equation, and gives back a `plot()` of the number of iteration vs. the recursive evaluations of the
logistic equation. `x0` is the initial condition, `r` is the rate parameter, n is the number of iterations, and `yn` is an optional argument,
when `true`, the function will plot the solutions.

"""

function logistic_iterator(x0, r, n, yn::Bool)

    solutions = []
    times = collect(0:n)

    push!(solutions, x0)

    for i in 1:n

        x_new = logistic_equation(x0, r)
        push!(solutions, x_new)
        x0 = x_new

    end

    plot(times, solutions)

    if yn == true

        plot(times, solutions, ".")

    end

    println(solutions[n+1])

end

"""
      iterator(f, n, x0)
      OUT:  steps, solutions

This is a function used to iterate a single variable function `f` (R -> R).

"""

function iterator(f::Function, n::Int, x0)

    solution = Float64[x0]
    steps = Int[0]
    x_old = x0
    for i in 1:n

        x_new = f(solution[i])
        push!(solution, x_new)
        push!(steps, i)

    end

    return steps, solution

end
