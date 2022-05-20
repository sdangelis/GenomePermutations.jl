using BED
using Distributions
using GenomicFeatures
using GenomePermutations
using HypothesisTests
using Statistics
using Test

@testset "GenomicPermutations.jl" begin
    
    # set up intervals to test
    interval1 = GenomicFeatures.Interval("chr1", 1150, 1180)  # test2
    interval2 = GenomicFeatures.Interval("chr1", 250,280) # test3
    interval3 = GenomicFeatures.Interval("chr1", 99,280) # test7
    interval4 = GenomicFeatures.Interval("chr3", 99,280) # test7
    interval5 = GenomicFeatures.Interval("chr2", 2150, 2180)
    interval6 = GenomicFeatures.Interval("chr2", 2350, 2380)
    collection1 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 100),
        GenomicFeatures.Interval("chr1", 200, 300),
        GenomicFeatures.Interval("chr1", 400, 500)
        ]) # test1
    collection2 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 50, 80),
        GenomicFeatures.Interval("chr1", 250, 280)
        ]) # test4
    collection3 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 50, 80),
        GenomicFeatures.Interval("chr1", 2250, 2280)
        ]) # test5 
    collection4 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 2150, 2180),
        GenomicFeatures.Interval("chr1", 2250, 2280)
        ]) # test6 
    collection5 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 2150, 2180),
        GenomicFeatures.Interval("chr1", 2250, 2280),
        GenomicFeatures.Interval("chr2", 2250, 2280)
        ])
    collection6 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr2", 2150, 2180),
        GenomicFeatures.Interval("chr2", 2250, 2280)
    ])
    collection7 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 0, 1000),
        GenomicFeatures.Interval("chr1", 2000, 3000)
    ])
    collection8 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr2", 2150, 2180),
        GenomicFeatures.Interval("chr2", 2250, 2280),
        interval6
    ])
    collection9 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 100),
        GenomicFeatures.Interval("chr1", 200, 300),
        GenomicFeatures.Interval("chr1", 2400, 2500),
    ])
    collection10 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 100),
        GenomicFeatures.Interval("chr1", 200, 300),
        GenomicFeatures.Interval("chr1", 2400, 2500),
        GenomicFeatures.Interval("chr2", 1, 100),
        GenomicFeatures.Interval("chr2", 200, 300),
        GenomicFeatures.Interval("chr2", 2400, 2500)
    ])
    collection11 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 5),
        GenomicFeatures.Interval("chr1", 21, 25),
        GenomicFeatures.Interval("chr1", 41, 45)
    ]) 
    collection12 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 5),
        GenomicFeatures.Interval("chr1", 21, 25),
        GenomicFeatures.Interval("chr1", 41, 45),
        GenomicFeatures.Interval("chr2", 1, 5),
        GenomicFeatures.Interval("chr2", 21, 25),
        GenomicFeatures.Interval("chr2", 41, 45),
    ]) 
    collection13 = GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1", 101, 199),
            GenomicFeatures.Interval("chr1", 501, 600),
            GenomicFeatures.Interval("chr1", 601, 701)
    ]) 
    collection14 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 100),
        GenomicFeatures.Interval("chr1", 200, 300),
        GenomicFeatures.Interval("chr1", 400, 500),
        GenomicFeatures.Interval("chr2", 1, 100),
        GenomicFeatures.Interval("chr2", 200, 300),
        GenomicFeatures.Interval("chr2", 400, 500)
    ])
    collection15 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 50, 80),
        GenomicFeatures.Interval("chr1", 250, 280),
        GenomicFeatures.Interval("chr2", 50, 80),
        GenomicFeatures.Interval("chr2", 250, 280)
    ]) 
    collection16 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 101, 199),
        GenomicFeatures.Interval("chr1", 501, 600),
        GenomicFeatures.Interval("chr1", 601, 701),
        GenomicFeatures.Interval("chr2", 101, 199),
        GenomicFeatures.Interval("chr2", 501, 600),
        GenomicFeatures.Interval("chr2", 601, 701)        
    ]) 
    collection17 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr2", 50, 80),
        GenomicFeatures.Interval("chr2", 250, 780)
    ]) 

    regions1 = collection7 
    regions2 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1", 0, 1000),
        GenomicFeatures.Interval("chr1", 2000, 3000),
        GenomicFeatures.Interval("chr2", 0, 1000),
        GenomicFeatures.Interval("chr2", 2000, 3000)
    ])
    regions3 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1",500,650),
        GenomicFeatures.Interval("chr1", 1000, 1101)
    ])
    regions4 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chr1",500,650),
        GenomicFeatures.Interval("chr1", 1000, 1101),
        GenomicFeatures.Interval("chr2",500,650),
        GenomicFeatures.Interval("chr2", 1000, 1101)
    ])
    regions5 = GenomicFeatures.IntervalCollection([
        GenomicFeatures.Interval("chrM", 0, 1000),
        GenomicFeatures.Interval("chrM", 2000, 3000),
    ])
    regions6 = IntervalCollection([GenomicFeatures.Interval("chr1", 1, 500000)])
    regions7 = IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 500000),
        GenomicFeatures.Interval("chr2", 1, 500000)
        ])
    regions8 = IntervalCollection([GenomicFeatures.Interval("chrM", 1, 500000)])  
    regions9 = IntervalCollection([GenomicFeatures.Interval("chr1", 1, 701)])
    regions10 = IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 701),
        GenomicFeatures.Interval("chr2", 1, 701)
    ])
    regions11 = IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 701),
        GenomicFeatures.Interval("chrM", 1, 701)
    ])
    regions12 =  IntervalCollection([
        GenomicFeatures.Interval("chr1", 1, 500000),
        GenomicFeatures.Interval("chr2", 1, 500000),
        GenomicFeatures.Interval("chrM", 1, 500000)
        ])

    # no need to test isoverlapping(::GenomicFeatures.Interval, ::GenomicFeatures.Interval)
    @testset "overlaps and isin functions" begin
        @test ([false, true, true, false] == map(x->anyoverlapping(x,collection1), [interval1, interval2, interval3, interval4]))
        @test GenomePermutations.alloverlapping(collection1, collection1)
        @test ([true, false, false]  == map(x->alloverlapping(x,collection1), [collection2, collection3, collection4]))
    
        @test isin(interval1,interval1) && isin(interval2,interval3) && !(isin(interval1,interval2)) && !(isin(interval3,interval4))
        @test isin(interval2, collection1) && !(isin(interval1, collection1))
        @test isin(collection1, collection1) && isin(collection2, collection1) && !(isin(collection3, collection2)) && !(isin(collection6, collection5))

        # NEED TO ADD A TEST FOR WHEN INTERVALS ARE contained in a collection --> need to differentiatie in/issubset and contain
        # test subsetting functions are equivalent 
        @test collection4 == GenomePermutations.vec_getcollection(collection5, "chr1") == GenomePermutations.iter_getcollection(collection5, "chr1")

        # test overlap counts
        @test countoverlapping(collection1, collection1) == collection1.length
        @test countoverlapping(collection2, collection1) == 2
        @test countoverlapping(collection6, collection1) == countoverlapping(collection1, collection6) == 0 
        @test countoverlapping(collection5, collection6) == 1
    end


    @testset "distance functions" begin
        # for single intervals
        @test ismissing(dist(interval3, interval4)) && dist(interval1, interval2) == 870 && dist(interval2, interval3) == 0
        
        # for an interval and a collection 
        @test [20,150,96.66666666666667,120] == map( x->dist(interval2 ,collection1, x), [minimum, maximum, mean, median])
        @test [970, 1070, 1020.0, 1020.0 ] == map( x->dist(interval1 ,collection5, x), [minimum, maximum, mean, median])
        #@test_logs (:warn, r"could not calculcate distance") match_mode=:any map( x->dist(interval1 ,collection5, x), [minimum, maximum, mean, median])
        
        # for 2 collections
        @test dist(collection5, collection5) == [0, 0, 0]
        @test  [[0, 170], 0, 170, 85.0, 85.0] == map( x->dist(collection2 ,collection3, x), [x -> x, minimum, maximum, mean, median])
        @test  [[2050, 1850, 1650], 1650, 2050, 1850.0, 1850.0] == map( x->dist(collection1 ,collection5, x), [x -> x, minimum, maximum, mean, median])
        #@test_logs (:warn, r"could not calculcate distance") match_mode=:any map( x->dist(collection1 ,collection5, x), [x -> x, minimum, maximum, mean, median])
        @test [[1650, 1750], 1650, 1750, 1700.0, 1700.0] ==  map( x->dist(collection5 ,collection1, x), [x -> x, minimum, maximum, mean, median]) # why does this fail
    end     


    # test distribution generation
    @testset "generate distributions" begin
        @test generatedistribution(collection1) == 
            MixtureModel([DiscreteUniform(1,100), DiscreteUniform(200, 300), DiscreteUniform(400, 500)], Distributions.Categorical([99/299, 100/299, 100/299]))
        @test typeof(generatedistribution(collection1)) <: Distributions.Distribution
    end 


    dist2 = generatedistribution(collection2)
    dist3 = generatedistribution(collection3)
    dist6 = generatedistribution(collection6)
    dist7 = generatedistribution(collection7)
    dist8 = generatedistribution(collection8)

    @testset "private interval randomisation" begin 
        @test_throws ErrorException GenomePermutations._randominterval(interval1, dist7, GenomicFeatures.IntervalCollection{nothing}(); allow_overlap = true)
        @test isin(GenomePermutations._randominterval(interval2, dist7, collection7; allow_overlap = true), collection7)
        @test isin(GenomePermutations._randominterval(interval2, dist7, collection7; allow_overlap = true), collection7) 
        @test_throws ErrorException GenomePermutations._randominterval(interval5, dist6, collection6; collection = collection6, allow_overlap = false)  # must throw exception
        @test isin(GenomePermutations._randominterval(interval5, dist6, collection8; collection = collection6, allow_overlap = true), collection8)
        @test GenomePermutations._randominterval(interval5, dist8, collection8; collection = collection6, allow_overlap = false) == interval6 # can only be 1 place
        @test_throws ErrorException GenomePermutations._randominterval(interval4, dist6, collection1; allow_overlap = false) # all stupidly mismatched. Indeed it fails gracefully
        @testset "fail behaviour" begin 
            # test default behaviour
            @test_throws ErrorException  GenomePermutations._randominterval(interval3, dist2, collection2; allow_overlap = false) 
            # test :orig
            @test GenomePermutations._randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :orig) == interval3
            # test :skip
            @test isnothing(GenomePermutations._randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :skip))
            # test :throw
            @test_throws ErrorException   GenomePermutations._randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :throw) 
        end 
    end 


    @testset "public interval randomisation" begin    
        # now we get prettier errors if sequences mismatch
        @test_throws KeyError GenomePermutations.randominterval(interval1, dist7, GenomicFeatures.IntervalCollection{nothing}(); allow_overlap = true)
        @test_nowarn GenomePermutations.randominterval(interval1, dist7, collection7; allow_overlap = true)
        @test isin(GenomePermutations.randominterval(interval2, dist7, collection7; allow_overlap = true, onfail = :orig), collection7)
        @testset "public fail behaviour" begin 
            # test default 
            @test_throws ErrorException  GenomePermutations.randominterval(interval3, dist2, collection2; allow_overlap = false) 
            # test :orig
            @test GenomePermutations.randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :orig) == interval3
            # test skip 
            @test isnothing(GenomePermutations.randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :skip))
            # test :throw
            @test_throws ErrorException   GenomePermutations.randominterval(interval3, dist2, collection2; allow_overlap = false, onfail = :throw) 
        end
    end


    @testset "region randomisation" begin 
        result = randomiseregions(collection1, dist7, collection7) 
        @test isin(result, collection7) && result != collection1 # check perfect case
        # check we fail if we must overlap
        @test_throws ErrorException randomiseregions(collection8, dist6, collection6)
        # check if allowing overlaps rescues test 2 
        @test_nowarn randomiseregions(collection8, dist6, collection6; allow_overlap = true)
        # check we converge to the only available solution
        @test randomiseregions(collection2, dist3, collection3) == collection3
        # check we converge to self
        @test randomiseregions(collection2, dist2, collection2; allow_overlap = false, max_tries=100000) == collection2
        #= This mess should erorr in an undefinend way. God only knows what kind of exception
        but it should do any error and refuse to go through an infinite loop =# 
        @test_throws ErrorException randomiseregions(collection3, dist2, collection6)
        @testset "check error handling" begin
            @test_throws ErrorException randomiseregions(collection17, dist6, collection6)
            @test_throws ErrorException randomiseregions(collection17, dist6, collection6, onfail = :throw)
            @test_throws ErrorException randomiseregions(collection17, dist6, collection6)
            @test randomiseregions(collection17, dist6, collection6; onfail = :skip).length == 1 
        end
    end


    @testset "genome randomisation" begin
        # test 1 sequence
        @test isin(randomisegenome(collection9, regions1), regions1)
        @test randomisegenome(collection9, regions1) != collection9
        # test 2 sequences
        @test isin(randomisegenome(collection10, regions2), regions2)
        @test randomisegenome(collection10, regions2) != collection10
        #= test over itself - note we need tiny interval to make sure we
        always draw the only possible solution within a reasonable number of tries =#
        @test randomisegenome(collection11, collection11; max_tries = 10000) == collection11
        @test randomisegenome(collection12, collection12; max_tries = 10000) == collection12
        # test a case that only works with no overlaps 
        @test_throws Core.Exception randomisegenome(collection9, regions3)
        @test isin(randomisegenome(collection9, regions3; allow_overlap=true), regions3)
        @test_throws Core.Exception randomisegenome(collection10, regions4)
        @test isin(randomisegenome(collection10, regions4; allow_overlap=true), regions4)
        # test an impossible case 
        @test_logs (:warn, r"is not found in the target regions")  randomisegenome(collection10, regions1) # skips wrong sequences 
        # all sequences are wrong - should return empty
        @test randomisegenome(collection9, regions5) == GenomicFeatures.IntervalCollection{Nothing}()
        
        @testset "onfail behaviour"  begin
            # test we successfully throw if we fail to randomise an interval
            @test_throws ErrorException randomisegenome(collection9, regions3; onfail = :throw) == collection9
            # test a case where we pass the original interval on fail
            @testset ":orig" begin
                scollection1 = GenomicFeatures.IntervalCollection([
                    GenomicFeatures.Interval("chr1", 1, 2),
                    GenomicFeatures.Interval("chr1",3, 4),
                    GenomicFeatures.Interval("chr1", 200, 500),
                    ]) 
                sregions1 = GenomicFeatures.IntervalCollection([
                    GenomicFeatures.Interval("chr1", 1, 4),
                    GenomicFeatures.Interval("chr1", 300, 400),
                    ]) 
                # test a case where we pass the original interval on fail and allow overlaps
                res = randomisegenome(scollection1, sregions1; onfail = :orig, allow_overlap = true)
                @test isin(GenomicFeatures.Interval("chr1", 200, 500), res)
                # really odd corner case
                @test_broken randomisegenome(scollection1, sregions1; onfail = :orig, allow_overlap = false) == scollection1
            end 
        end
    end
    
     
    @testset "simple P structure (constructor)" begin 
        sp = SimplePTest(100, 0.1, :Less) 
        @test sp.iterations == 100
        @test sp.p_val == 0.1
        @test sp.alternative == :Less
        @test sp.min_p_val == 0.009900990099009901
        @test pvalue(sp) == 0.1
    end 


    @testset "PermTestResult structure (constructor)" begin 
        @testset "single value" begin 
            res = PermTestResult(3, 10, [10, 10, 10], SimplePTest(3, 1, :Less) , "test1", "test2", overlappermtest ,"guess and check")
            @test typeof(res) == PermTestResult{Int, Int, SimplePTest}
            @test length(res.ran) == res.iterations
            @test pvalue(res) == 1.0
        end 
        @testset "array" begin 
            res = GenomePermutations.PermTestResult(3, [1, 2, 3], [[1, 2, 3],[1, 2, 3],[1, 2, 3]], SimplePTest(3, 1, :Less) , "test1", "test2", overlappermtest ,"guess and check")
            @test typeof(res) == PermTestResult{Vector{Int}, Vector{Int}, SimplePTest}
            @test length(res.ran) == res.iterations
            @test pvalue(res) == 1.0
        end
    end 
    
    
    @testset "simple P test" begin
        @test GenomePermutations.simpleP(80, [80, 80, 80, 80, 80], 5; alternative = "Auto") == GenomePermutations.SimplePTest(5, 1, :Less)
        @test GenomePermutations.simpleP(80, [400, 400, 400, 400, 400, 400, 400, 400, 400, 400], 10; alternative = "Auto") == GenomePermutations.SimplePTest(10, 0.09090909090909091, :Less)
        @test GenomePermutations.simpleP(80, [78, 78, 78, 99, 78], 5; alternative = "Auto") == GenomePermutations.simpleP(80, [2, 2, 2, 99, 2], 5; alternative = "Auto") ==GenomePermutations.SimplePTest(5, 0.3333333333333333,:More)
    end 

    
    @testset "overlap permutation test" begin 
        @testset "single sequence" begin
            # 1 Overlap with itself should be as significant as it can
            @test overlappermtest(collection1,  collection1, regions6, 100).test.p_val == 
                GenomePermutations.simpleP(1,repeat([1], 100),100).min_p_val
            # 2 Less overlaps than expected
            test = overlappermtest(collection1, collection13, regions9, 100)
            @test test.test.p_val < 0.01 && test.test.alternative == :Less
            # 3 More overlaps than expected
            test = overlappermtest(collection1, collection2, regions6, 100)
            @test test.test.p_val < 0.01 && test.test.alternative == :More
            # 4 sequenecs to randomise do not belong to regions 
            @test_throws ErrorException overlappermtest(collection1,  collection2, regions8, 100)
            # 5 Regions to randomise have an extra sequence
            @test_nowarn overlappermtest(collection1,  collection2, regions11, 100)            
        end
        @testset "multiple sequences" begin 
            # 1 Overlap with itself should be as significant as it can
            @test overlappermtest(collection10,  collection10, regions7, 100).test.p_val == 
                GenomePermutations.simpleP(1,repeat([1], 100),100).min_p_val
            # 2 Less overlaps than expected
            test = overlappermtest(collection14,  collection16, regions10, 100)
            @test test.test.p_val < 0.01 && test.test.alternative == :Less
            # 3 More overlaps than expected
            test = overlappermtest(collection14,  collection15, regions7, 100)
            @test test.test.p_val < 0.01 && test.test.alternative == :More
            # 4 All sequences wrong
            @test_throws ErrorException overlappermtest(collection14,  collection15, regions8, 100)
            # 5 One of 2 sequences are wrong
            @test_throws ErrorException overlappermtest(collection14,  collection15, regions11, 100)
            # 6 there is an extra chromosome - This should work fine
            @test_nowarn overlappermtest(collection14,  collection15, regions12, 100)
        end
    end
    
    @testset "Random Genome Generation" begin
        @test_nowarn GenomePermutations._randomgenome(regions2, 10, Distributions.Categorical(2), Distributions.Binomial(10),["chr1","chr2"])
        @test isin(GenomePermutations._randomgenome(regions2, 10, Distributions.Categorical(2), Distributions.Binomial(10),["chr1","chr2"]), regions2)
        @test_throws KeyError GenomePermutations._randomgenome(regions1, 10, Distributions.Categorical(2), Distributions.Binomial(10),["chr1","chr2","chrZ"])
        @test_throws KeyError GenomePermutations._randomgenome(regions1, 10, Distributions.Categorical(1), Distributions.Binomial(10),["chr1","chr2","chrZ"])
    end
    
    @testset "Generic permutation test" begin 
        # define some supplementary regions as needed 
        scollection1 =  GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1", 2, 101), 
            GenomicFeatures.Interval("chr1", 201, 301), 
            GenomicFeatures.Interval("chr1", 402, 503)])
        scollection2 =  GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1", 2, 101), 
            GenomicFeatures.Interval("chr1", 201, 301), 
            GenomicFeatures.Interval("chr1", 402, 503),
            GenomicFeatures.Interval("chr2", 2, 101), 
            GenomicFeatures.Interval("chr2", 201, 301), 
            GenomicFeatures.Interval("chr2", 402, 503)
            ])
        scollection3 =  GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1",  499400, 499500), 
            GenomicFeatures.Interval("chr1",  499700, 499800), 
            GenomicFeatures.Interval("chr1", 499900, 500000)])
        scollection4 =  GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1",  499400, 499500), 
            GenomicFeatures.Interval("chr1",  499700, 499800), 
            GenomicFeatures.Interval("chr1", 499900, 500000),
            GenomicFeatures.Interval("chr2",  499400, 499500), 
            GenomicFeatures.Interval("chr2",  499700, 499800), 
            GenomicFeatures.Interval("chr2", 499900, 500000)])
        scollection5 = GenomicFeatures.IntervalCollection([
                GenomicFeatures.Interval("chr1", 1, 100),
                GenomicFeatures.Interval("chr1", 200, 300),
                GenomicFeatures.Interval("chr1", 400, 500),
                GenomicFeatures.Interval("chr2", 1, 100),
                GenomicFeatures.Interval("chr2", 200, 300),
                GenomicFeatures.Interval("chr2", 400, 500)
                ]) 

        @testset "default distance" begin
            @testset "Single sequence" begin 
                # 1 distance to iteslf should be as significant as it gets
                test =  permtest(collection1, collection1, regions6, 100) 
                @test pvalue(test.test) < 0.01
                # 2 closer than expected
                test = permtest(collection1, scollection1, regions6, 100) 
                @test pvalue(test.test) < 0.01
                # 3 more distant than expected
                test = permtest(collection1, scollection3, regions6, 100)
                @test pvalue(test.test) < 0.01
                # 4 sequenecs to randomise do not belong to regions 
                @test_throws ErrorException permtest(collection1,  collection2, regions8, 100)
                # 5 Regions to randomise have an extra sequence
                @test_nowarn permtest(collection1,  collection2, regions11, 100)  
            end
            @testset "Mutliple sequences" begin 
                # 1 distance to iteslf should be as significant as it gets
                test =  permtest(collection10, collection10, regions7, 100) 
                @test pvalue(test.test) < 0.01
                # 2 closer than expected
                test = permtest(collection10, scollection2, regions7, 100) 
                @test pvalue(test.test) < 0.01
               # 3 more distant than expected
                test = permtest(scollection5, scollection4, regions7, 100)
                @test pvalue(test.test) < 0.01
                # Need a test for the case where we expect a p of 1. 
                # 4 All sequences wrong
                @test_throws ErrorException permtest(collection14,  collection15, regions8, 100)
                # 5 One of 2 sequences are wrong
                @test_throws ErrorException permtest(collection14,  collection15, regions11, 100)
                # 6 there is an extra chromosome - This should work fine
                @test_nowarn permtest(collection14,  collection15, regions12, 100)
            end

        end

        @testset "Large Collection test" begin 
            lr1 = IntervalCollection(sort!([Interval(i, 1, 100000) for i in "chr" .* map(string, collect(1:22))]))
            lr2 = IntervalCollection(sort!([Interval(i, 1, 100000) for i in "chr" .* map(string, collect(1:22))]))
            lc1 = GenomePermutations._randomgenome(lr1, 400, Distributions.Categorical(22), Distributions.Binomial(500),"chr" .* map(string, collect(1:22)))
            lc2 = GenomePermutations._randomgenome(lr2, 400, Distributions.Categorical(22), Distributions.Binomial(500),"chr" .* map(string, collect(1:22)))
            @test isin(lc1, lr1) &&  isin(lc2, lr2)
            @test_nowarn permtest(lc1,  lc2, lr1, 1000)
        end
    end

end
