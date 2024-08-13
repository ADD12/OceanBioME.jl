#Checked forcing equation, still to check parameters correctly named
#This document contains functions for:
    #λₚₛᵢ¹ (eq52, parametrisation of dissolution rate of PSi) 
    #Forcing for PSi (eq51)

@inline function χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc)

    χ_lab⁰ = bgc.proportion_of_the_most_labile_phase_in_PSi
    zₘₐₓ = max(abs(zₘₓₗ), abs(zₑᵤ))

    if abs(z) <= zₘₐₓ
        return χ_lab⁰
    else
        return χ_lab⁰*exp(-(λₚₛᵢˡᵃᵇ - λₚₛᵢʳᵉᶠ)*((abs(z)-zₘₐₓ)/(get_w_GOC(z, zₑᵤ, zₘₓₗ, bgc) + eps(0.0)))) #eq53
    end
end

@inline function get_λₚₛᵢ¹(zₘₓₗ, zₑᵤ, z, T, Si, bgc)

    λₚₛᵢˡᵃᵇ = bgc.fast_dissolution_rate_of_PSi
    λₚₛᵢʳᵉᶠ = bgc.slow_dissolution_rate_of_PSi

    Si_eq = 10^(6.44 - 968/(T + 273.15)) 
    Siₛₐₜ = (Si_eq - Si)/(Si_eq + eps(0.0))

    λₚₛᵢ = χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc)*λₚₛᵢˡᵃᵇ + (1 - χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc))*λₚₛᵢʳᵉᶠ

    return λₚₛᵢ*(0.225*(1 + T/15)*Siₛₐₜ + 0.775*(((1 + T/400)^4)*Siₛₐₜ)^9) 
end

@inline function (bgc::PISCES)(::Val{:PSi}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³) 
    #Parameters
    Kₘ = bgc.half_saturation_const_for_mortality
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    αᴰ = bgc.initial_slope_of_PI_curve.D
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    Dissₛᵢ = bgc.dissolution_rate_of_silicon
    mᴰ = bgc.phytoplankton_mortality_rate.D
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton

    #Also required
    sh = get_sh(z, zₘₓₗ)
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ)
    θˢⁱᴰ = θ(Dˢⁱ, D)

    return  (θˢⁱᴰ*get_grazingᴹ(P, D, Z, POC, T, bgc)[3]*M +  θˢⁱᴰ*get_grazingᶻ(P, D, POC, T, bgc)[3]*Z 
           + mᴰ*concentration_limitation(D, Kₘ)*Dˢⁱ + sh*wᴰ*D*Dˢⁱ - get_λₚₛᵢ¹(zₘₓₗ, zₑᵤ, z, T, Si, bgc)*Dissₛᵢ*PSi) #removed θˢⁱᴰ from third term, to conserve silicon
end