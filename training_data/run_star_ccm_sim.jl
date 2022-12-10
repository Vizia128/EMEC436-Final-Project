include("dat_to_sdf.jl")


function modify_java_macro(name, naca, x, y, magnitude)
    open("record_variable1.java", "r+") do file
        lines = readlines(file)

        for (i, line) in enumerate(lines)
            if i == 58
                lines[i] = """    cadModel_0.getFeatureManager().create3DSketches_2("C:\\\\Users\\\\kizan\\\\OneDrive\\\\Documents\\\\__MSU\\\\Fall 2022\\\\EGEN 436\\\\Final Project\\\\Airfoils\\\\CSV_naca\\\\NACA$naca.csv", labCoordinateSystem_0, false, true);"""
            end
            if i == 195
                lines[i] = """    flowDirectionProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits($x, $y, 0.0, units_2);"""
            end
            if i == 203
                lines[i] = """    velocityMagnitudeProfile_0.getMethod(ConstantScalarProfileMethod.class).getQuantity().setValueAndUnits($magnitude, units_3);"""
            end
            if i == 211
                lines[i] = """    velocityProfile_0.getMethod(ConstantVectorProfileMethod.class).getQuantity().setComponentsAndUnits($(x*magnitude), $(y*magnitude), 0.0, units_3);"""
            end
            if i == 231
                lines[i] = """    tableUpdate_0.setBaseFilename("$name, all");"""
            end
            if i == 239
                lines[i] = """    tableUpdate_1.setBaseFilename("$name, rect");"""
            end
            if i == 254
                lines[i] = """    plotUpdate_0.setAnimationFilenameBase("$name, res");"""
            end
        end

        println(lines[58])

        write("star_macro_run.java", join(lines, "\n"))
    end
end

function get_random_parameters()
    AoA = rand()
    AoA = (AoA - 0.5) * 32

    logRe = rand()
    logRe = logRe*4 + 2 
    Re = exp10(logRe)

    n1 = rand(0:6)
    n2 = rand(0:6)
    n34 = rand(6:24)
    n34 = lpad(n34, 2, '0')

    naca = "$n1$n2$n34"

    return AoA, Re, naca
end

function run_simulation(;n=1)
    for i in 1:n
        # cd("C:/Users/kizan/OneDrive/Documents/__MSU/Fall 2022/EGEN 436/Final Project")
        AoA, Re, naca = get_random_parameters()

        x = cos(AoA*pi/180)
        y = sin(AoA*pi/180)

        μ = 1.85508e-5 # Pa-s
        ρ = 1.18415 # kg/m^3
        velocity_magnitude = μ*Re/ρ

        modify_java_macro("naca$naca, AoA=$AoA, Re=$Re", naca, x, y, velocity_magnitude)

        Airfoil("Airfoils/naca/naca$(naca).dat")

        # cd("C:/Users/kizan/OneDrive/Documents/__MSU/Fall 2022/EGEN 436/Final Project/StarCCM")
        run(`starccm+ -np 12 -power -podkey 2PBihIVQfKVNEJmuUmShHw -licpath 1999@flex.cd-adapco.com -batch star_macro_run.java -load airfoil_flow.sim`)
    end
end
@time run_simulation(n=300)