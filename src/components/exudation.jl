function exudation_io(filestr)
    # filestr = /files/TOC_data_Hopland_Avena.csv
    dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
    df = CSV.read(filestr, DataFrame, missingstring = "N/A")
    root_density  = 0.6e-3                                                      # g/g, Rillig (2002)
    bulk_density  = 1.21                                                        # g/cm^3, Thea Whitman
    week3tot      = median(df[1:4, 2])*root_density*bulk_density
    week6tot      = median(df[5:8,2])*root_density*bulk_density
    week9tot      = median(df[9:12,2])*root_density*bulk_density
    week12tot     = median(df[13:end,2])*root_density*bulk_density
    TOC           = [0.0, week3tot, week6tot, week9tot, week12tot, 0.0].*12.011      # μg/cm^3 -> mol/m^3
    x0            = [0, 3*7.0, 6*7.0, 9*7.0, 12*7.0, 15*7.0].*24                    # hours
    return x0, TOC
end

function exudation_class_io(filestr)
    dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
    df = CSV.read(joinpath(dir, filestr), DataFrame, missingstring = "N/A")
    df = filter(row -> row[:Ontology] !== "other", df)
    df_class  = @linq df |>
                      by(:Ontology, week3 = mean(:week3_rnorm), week6 = mean(:week6_rnorm),
                                    week9 = mean(:week9_rnorm), week12 = mean(:week12_rnorm)) |>
                      transform(week3_norm = :week3./sum(:week3), week6_norm = :week6./sum(:week6),
                                week9_norm = :week9./sum(:week9), week12_norm = :week12./sum(:week12))


    class_norm = convert(Matrix, df_class[:,[:week3_norm, :week6_norm, :week9_norm, :week12_norm]])
    amino_class = class_norm[1,:]
    organic_class = class_norm[2,:]
    nucleotide_class = class_norm[3,:]
    sugar_class = class_norm[4,:]
    auxin_class = class_norm[5,:]
    fatty_class = class_norm[6,:]
    return  hcat(zeros(6), class_norm, zeros(6))
end

function max_exudation_rate(xt, TOC, t)
    spl = Dierckx.Spline1D(xt, TOC)
    year_t = 0:365*24
    tmp = max.(0.0, evaluate(spl, year_t))
    spl2 = Spline1D(year_t[1:end-1], diff(tmp)./diff(year_t), s=1e-11)
    return evaluate(spl2, t)
end
