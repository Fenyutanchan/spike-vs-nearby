# Copyright (c) 2026 Quan-feng Wu <wqf@fytc.ac>
#
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

using QuadGK: quadgk

@testset "homogeneous transport physics" begin
    D0 = 3.0e28 * NU.cm^2 / NU.s
    b0 = 1.0e-16 * EU() / NU.s
    diffusion = DiffusionModel_PowerLaw(D0, 1 / 3, EU(1.0))
    losses = EnergyLossModel_PowerLaw(b0, 2.0, EU(1.0))
    transport = TransportModel_Homogeneous(diffusion, losses)

    @test diffusion_coefficient(transport, EU(1.0)) / D0 ≈ 1.0
    @test diffusion_coefficient(transport, EU(1.0e3)) / D0 ≈ 10.0
    @test energy_loss_rate(transport, EU(1.0)) / b0 ≈ 1.0
    @test energy_loss_rate(transport, EU(10.0)) / b0 ≈ 100.0

    E1, E2, E3 = EU(10.0), EU(100.0), EU(1.0e4)
    numerical_time, _ = quadgk(10.0, 1.0e4) do energy
        EUval(EU, inv(energy_loss_rate(transport, EU(energy))))
    end
    numerical_integral, _ = quadgk(10.0, 1.0e4) do energy
        EUval(
            EU,
            diffusion_coefficient(transport, EU(energy)) /
            energy_loss_rate(transport, EU(energy)),
        )
    end
    @test unit_isapprox(
        cooling_time(transport, E3, E1),
        EU(numerical_time, -1);
        rtol=1e-12,
    )
    @test unit_isapprox(
        diffusion_loss_integral(transport, E3, E1),
        EU(numerical_integral, -2);
        rtol=1e-12,
    )
    @test unit_isapprox(
        maximum_cooling_time(transport, E1),
        E1 / energy_loss_rate(transport, E1),
    )

    elapsed_time = cooling_time(transport, E3, E2)
    @test unit_isapprox(
        source_energy_ceiling(transport, E2, elapsed_time),
        E3;
        rtol=1e-12,
    )
    @test isinf(source_energy_ceiling(
        transport,
        E2,
        maximum_cooling_time(transport, E2),
    ))

    @test unit_isapprox(
        maximum_diffusion_loss_integral(transport, E1),
        E1 * diffusion_coefficient(transport, E1) /
        energy_loss_rate(transport, E1) / (2 - 1 / 3 - 1),
    )

    diffusion_integral = diffusion_loss_integral(transport, E3, E2)
    upper_radius = 15 * sqrt(EUval(EU, diffusion_integral))
    kernel_integral, _ = quadgk(0.0, upper_radius) do radius
        4 * pi * radius^2 *
        EUval(EU, green_kernel(diffusion_integral, EU(radius, -1)))
    end
    second_moment, _ = quadgk(0.0, upper_radius) do radius
        4 * pi * radius^4 *
        EUval(EU, green_kernel(diffusion_integral, EU(radius, -1)))
    end
    @test kernel_integral ≈ 1.0 rtol=1e-12
    @test second_moment ≈ 6 * EUval(EU, diffusion_integral) rtol=1e-12
end
