@testset "quasi-static-inclined-support" begin
    cp("../examples/single/static-solid/cube_compression_inclined_support_30_degrees/cube.yaml", "cube.yaml", force=true)
    cp("../examples/single/static-solid/cube_compression_inclined_support_30_degrees/cube.g", "cube.g", force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    
    # These displacements are obtained from an identical problem with
    # BCs coincident with the lab basis
    reference_displacements = [
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.19999999999999996, 0.0, 0.0, 
        -0.002490510706379005, -0.054808932034627424, 
        0.054808932034627285, 0.008794796508373306, 
        -0.05449689589044056, -2.8756551738792074e-17, 
        0.020548157939405773, -7.567237356520244e-17,
        -4.131059890383733e-18, 0.008794796508373347, 
        -4.2461283759579046e-17, 0.05449689589044062, 
        0.19999999999999996, 0.0, 0.0, 
        0.19999999999999996, 0.0, 0.0, 
        -0.0024905107063790998, -0.05480893203462743, 
        -0.05480893203462725, 0.008794796508373271,
        -4.832171592315198e-17, -0.054496895890440634, 
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 
        0.0, 0.0, 0.008794796508373334, 
        0.05449689589044063, -5.145650465803633e-17,
        -0.002490510706379087, 0.05480893203462731, 
        0.0548089320346273, 0.19999999999999996, 0.0,
        0.0, -0.002490510706378977, 
        0.05480893203462727, -0.05480893203462724, -0.2,
        -0.052612760978146725, 
        0.052612760978146586, -0.2, -0.05432895977117359,
        9.641660437651523e-18, -0.2, 
        -3.8460049041556104e-17, -3.566854429999362e-17, 
        -0.2, -7.231031457420035e-17, 
        0.05432895977117351, -0.2, -0.05261276097814662,
        -0.052612760978146676, -0.2, 
        -1.0483851671759644e-16, -0.05432895977117339,
        -0.2, 0.054328959771173434, 
        -2.6230761947053614e-17, -0.2, 0.052612760978146704,
        0.05261276097814651, -0.2, 
        0.05261276097814652, -0.05261276097814667
    ]

    # Rotate these displacements
    angle = 30. * π / 180
    c = cos(angle)
    s = sin(angle)
    # Rotate about z
    local_rotation_matrix = [ c s 0; -s c 0 ; 0 0 1]
    global_rotation = zeros((81,81))
    for i in range(1,27)
        base = (i-1)*(3) + 1
        global_rotation[base:(base+2), base:(base+2)] = local_rotation_matrix
    end

    correct_displacements = global_rotation' * reference_displacements

    local_stiffness_displacements = integrator.displacement
    # Displacements in the integrator are stored in the local frame
    # (i.e., inclined support frame). Use the model global_transform
    # to return them to the global frame
    displacements = model.global_transform' * integrator.displacement

    # Assert the displacement array matches the reference displacements
    @test displacements ≈ correct_displacements atol=1e-6
    
end

@testset "newark-inclined-support" begin
    cp("../examples/single/implicit-dynamic-solid/cube_compression_inclined_support_30_degrees/cube.yaml", "cube.yaml", force=true)
    cp("../examples/single/implicit-dynamic-solid/cube_compression_inclined_support_30_degrees/cube.g", "cube.g", force=true)
    simulation = Norma.run("cube.yaml")
    integrator = simulation.integrator
    model = simulation.model
    rm("cube.yaml")
    rm("cube.g")
    rm("cube.e")
    
    # These displacements are obtained from an identical problem with
    # BCs coincident with the lab basis
    reference_displacements = reference_displacements = [
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 0.0, 0.0,
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 0.0, 0.0,
        -0.002490463082911527, -0.05480993185221903, 0.05480993185221887,
        0.008795258511939654, -0.05449788272457724, -3.894904718685741e-17,
        0.02054905860583953, -5.1586464370647656e-17, -1.3181110813082193e-17,
        0.008795258511939637, -1.7537207823376988e-17, 0.054497882724577194,
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 0.0, 0.0,
        -0.0024904630829116706, -0.05480993185221884, -0.05480993185221892,
        0.008795258511939583, 1.6885023924638925e-17, -0.05449788272457717,
        0.19999999999999996, 0.0, 0.0, 0.19999999999999996, 0.0, 0.0,
        0.008795258511939656, 0.05449788272457725, 1.2228380319626057e-17,
        -0.002490463082911655, 0.05480993185221885, 0.054809931852218774,
        0.19999999999999996, 0.0, 0.0, -0.0024904630829117136, 0.05480993185221884,
        -0.054809931852218774, -0.2, -0.052614066999186214, 0.0526140669991862,
        -0.2, -0.05433034200880138, -6.531259792682969e-17, -0.2,
        2.133495509760114e-17, -7.3511277735106e-17, -0.2, -1.167349277228399e-16,
        0.054330342008801345, -0.2, -0.05261406699918625, -0.05261406699918614,
        -0.2, -3.366299465092707e-17, -0.054330342008801394, -0.2,
        0.05433034200880128, 1.0918779448161045e-17, -0.2, 0.05261406699918612,
        0.0526140669991861, -0.2, 0.052614066999186256, -0.052614066999186124
    ]
    
    # Rotate these displacements
    angle = 30. * π / 180
    c = cos(angle)
    s = sin(angle)
    # Rotate about z
    local_rotation_matrix = [ c s 0; -s c 0 ; 0 0 1]
    global_rotation = zeros((81,81))
    for i in range(1,27)
        base = (i-1)*(3) + 1
        global_rotation[base:(base+2), base:(base+2)] = local_rotation_matrix
    end

    correct_displacements = global_rotation' * reference_displacements

    local_stiffness_displacements = integrator.displacement
    # Displacements in the integrator are stored in the local frame
    # (i.e., inclined support frame). Use the model global_transform
    # to return them to the global frame
    displacements = model.global_transform' * integrator.displacement

    # Assert the displacement array matches the reference displacements
    @test displacements ≈ correct_displacements atol=1e-5
    
end