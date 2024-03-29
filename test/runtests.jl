using BED
using Distributions
using GenomicFeatures
using GenomePermutations
using HypothesisTests
using Statistics
using Test

# set up intervals to test

interval1 = GenomicFeatures.Interval("chr1", 1150, 1180)  # test2
interval2 = GenomicFeatures.Interval("chr1", 250,280) # test3
interval3 = GenomicFeatures.Interval("chr1", 99,280) # test7
interval4 = GenomicFeatures.Interval("chr3", 99,280) # test7
interval5 = GenomicFeatures.Interval("chr2", 2150, 2180)
interval6 = GenomicFeatures.Interval("chr2", 2350, 2380)
interval8 = GenomicFeatures.Interval("chr2", 1350, 2400)
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

collection18 = GenomicFeatures.IntervalCollection([
    GenomicFeatures.Interval("chr2", 1, 5),
    GenomicFeatures.Interval("chr2", 21, 25),
    GenomicFeatures.Interval("chr2", 41, 45),
    GenomicFeatures.Interval("chr2", 101, 105),
    GenomicFeatures.Interval("chr2", 121, 125),
    GenomicFeatures.Interval("chr2", 141, 145),
    GenomicFeatures.Interval("chr2", 201, 205),
    GenomicFeatures.Interval("chr2", 221, 225),
    GenomicFeatures.Interval("chr2", 241, 245),
    GenomicFeatures.Interval("chr2", 301, 305),
    ]) 
collection20 = GenomicFeatures.IntervalCollection([
    GenomicFeatures.Interval("chr1", 101, 199),
    GenomicFeatures.Interval("chr1", 121, 199),
    GenomicFeatures.Interval("chr1", 501, 600),
    GenomicFeatures.Interval("chr1", 521, 600),
    GenomicFeatures.Interval("chr1", 601, 701),    
    GenomicFeatures.Interval("chr1", 621, 701)        
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
regions13 = IntervalCollection([
    GenomicFeatures.Interval("chr2", 1, 500000),
    GenomicFeatures.Interval("chr3", 1, 500000),
    GenomicFeatures.Interval("chrM", 1, 500000)
    ])
regions14 = IntervalCollection([
    GenomicFeatures.Interval("chr1", 1, 180),
    GenomicFeatures.Interval("chr1", 400, 550),
    GenomicFeatures.Interval("chr1", 700, 1000)
    ])


@testset "GenomicPermutations.jl" begin
    
    @testset "exports" begin 
    end 

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
        @test countoverlaps(collection1, collection1) == length(collection1)
        @test countoverlaps(collection20, collection20) > length(collection20)
        @test countoverlaps(collection20, collection20) > countoverlapping(collection20, collection20)
        @test countoverlapping(collection2, collection1) == 2
        @test countoverlapping(collection6, collection1) == countoverlapping(collection1, collection6) == 0 
        @test countoverlapping(collection5, collection6) == 1
    end

    @testset "distance functions" begin
        @testset "simple distance" begin 
            @testset "single interval"  begin
                @test ismissing(dist(interval3, interval4)) && dist(interval1, interval2) == 870 && dist(interval2, interval3) == 0
            end 
            @testset "single interval and collection" begin 
                @test 20 == dist(interval2 ,collection1)
                @test 970 == dist(interval1 ,collection5)
            end
            @testset "for 2 collections" begin
                @test dist(collection5, collection5) == [0, 0, 0]
                @test  [0, 170] == dist(collection2 ,collection3)
                @test  [2050, 1850, 1650] == dist(collection1 ,collection5)
                @test [1650, 1750] ==  dist(collection5 ,collection1) 

            end
        end

        @testset "feature distance" begin
            # if no overlaps featuredist == dist, else = 0
            @testset "single interval" begin 
                @test ismissing(featuredist(interval3, interval4)) && featuredist(interval1, interval2) == 870 && featuredist(interval2, interval3) == 0
                @test featuredist(interval6, interval8) == 0 != dist(interval6, interval8)
            end
            @testset "single interval and a feature collection" begin 
                @test 970 == dist(interval1 ,collection5) == featuredist(interval1 ,collection5)
                @test 0 == featuredist(interval2 ,collection1) != dist(interval2 ,collection1)
            end 
            @testset "interval and feature collection" begin 
                @test featuredist(collection5, collection5) == [0, 0, 0]
                @test featuredist(collection2, collection1) == [0,0] != dist(collection2, collection1)
                # test a partially impossible case
                @test featuredist(collection12, collection9) == [0, 0, 0] != dist(collection12, collection9)
                @test featuredist(collection12, regions14) == [0, 0, 0] != dist(collection16, regions14)
                @test featuredist(regions4, regions2) == [0, 0, 0, 0] != dist(regions4, regions2)
            end 
            # case where we have mno verla 
            #case where we have overlaps
        end
    end     

    @testset "test AbstractGenomeRand" begin
        
        @testset "StartMixture" begin
            
            @test_nowarn  GenomePermutations.StartMixture("hgTest", collection1, false) 
            @test_nowarn  GenomePermutations.StartMixture("hgTest", collection1, true) 
            
            t1 = GenomePermutations.StartMixture("hgTest", collection1, true) 
            t7 = GenomePermutations.StartMixture("hgTest", collection7, true)
            t17 = GenomePermutations.StartMixture("hgTest", collection17, false) 
            t17_by_chr = GenomePermutations.StartMixture("hgTest", collection17, true; by_chromosome = true) #
            t6 = GenomePermutations.StartMixture("hgTest", collection6, true; max_tries = 1000)
            f6 = GenomePermutations.StartMixture("hgTest", collection6, false)
            f8 = GenomePermutations.StartMixture("hgTest", collection8, false; max_tries = 100000)
            t8 = GenomePermutations.StartMixture("hgTest", collection8, true; max_tries = 100000)
            f8_by_chr = GenomePermutations.StartMixture("hgTest", collection8, false; by_chromosome = true, max_tries = 100000)
            f1 = GenomePermutations.StartMixture("hgTest", collection1, false) 
            t2r = GenomePermutations.StartMixture("hgTest", regions2, true) 
            f2r = GenomePermutations.StartMixture("hgTest", regions2, false)
            @test t1._distribution["chr1"] == MixtureModel([DiscreteUniform(1,100), DiscreteUniform(200, 300), DiscreteUniform(400, 500)], Distributions.Categorical([99/299, 100/299, 100/299]))
            
            @testset "draw random samples" begin
                @test_nowarn rand(f1)
                @test_nowarn rand(t1)
                @test isin(rand(f1), collection1)
                @test isin(rand(t1), collection1)
            end    
            
            @testset "randomiser" begin
                @test_nowarn  GenomePermutations.randomise(interval1, t17) 
                @test_throws KeyError GenomePermutations.randomise(interval1, t17_by_chr)
                @test isin(GenomePermutations.randomise(interval1, t7), collection7)
                # no way to randomise this interval without overlapping collection 6 
                @test_throws ErrorException GenomePermutations.randomise(interval5, f6; collection = collection6)        
                f6orig = deepcopy(f6)
                f6orig.on_fail = :orig
                @test GenomePermutations.randomise(interval5, f6orig; collection = collection6) == interval5
                @test isin(GenomePermutations.randomise(interval5, t6; collection = collection6), collection6)
                # this one can only be in 1 place 
                @test GenomePermutations.randomise(interval5, f8; collection = collection6) == interval6
                # are chromosome by chromosome the same when they must be?
                GenomePermutations.randomise(interval5, f8; collection = collection6) == \
                GenomePermutations.randomise(interval5, f8_by_chr; collection = collection6)
            end 

            @testset "Abstract Methods" begin
                # Do we work with the abstract methods?
                @testset "rand" begin
                @test_nowarn rand(t1,100) 
                @test typeof(rand(f1, 1)) <: GenomicFeatures.IntervalCollection
                # i with a small choice of locations we can't draw lots without overlapping
                @test_nowarn rand(f1, 1)
                @test_throws ErrorException rand(f1, 100)
                end

                @testset "randomise" begin            
                t13_by_chr = GenomePermutations.StartMixture("hgTest", regions13, true; by_chromosome = true)
                t13 = GenomePermutations.StartMixture("hgTest", regions13, true; by_chromosome = false)
                
                # inoffensive case 
                @test isin(GenomePermutations.randomise(collection11, t1), collection1)
                @test isin(GenomePermutations.randomise(collection11, f1), collection1)
                
                # check we don't overlap unless we ask to
                
                s1 = GenomicFeatures.IntervalCollection([
                    GenomicFeatures.Interval("chr1", 1, 50),
                    GenomicFeatures.Interval("chr1", 101, 150),
                    GenomicFeatures.Interval("chr1", 201, 250),
                    GenomicFeatures.Interval("chr1", 350, 400), 
                    GenomicFeatures.Interval("chr1", 400, 450),
                    GenomicFeatures.Interval("chr1", 500, 550)
                ])
                sc1 = GenomicFeatures.IntervalCollection([GenomicFeatures.Interval("chr1", 1, 500)])
                tsc1 = GenomePermutations.StartMixture("hgTest", sc1, true)
                fsc1 = GenomePermutations.StartMixture("hgTest", sc1, false, max_tries = 1000000)

                overlaps = GenomePermutations.randomise(s1, tsc1)
                no_overlaps = GenomePermutations.randomise(s1, fsc1)
                @test isin(overlaps, regions2)
                @test isin(no_overlaps, regions2) # we're in the target regions
                @test GenomePermutations.countoverlaps(overlaps, overlaps) > length(overlaps)  # if we allow overlaps more will overlap 
                @test GenomePermutations.countoverlaps(no_overlaps, no_overlaps) == length(no_overlaps) # if we disallow overlaps we only overlap with ourselves
                # test by chromosome returns always something on the right chromosome 
                @test collect(keys(GenomePermutations.randomise(collection18, t13_by_chr).trees)) == ["chr2"]
                # sequence mismatch -> does only work wihtout by_chromosome
                @test_nowarn GenomePermutations.randomise(collection5, t13)
                @test_throws KeyError GenomePermutations.randomise(collection5, t13_by_chr)
                # can only work without overlaps
                c6c6= deepcopy(collection6) # we don't know what uses collection 6 later and don't wanna mess it
                for interval in collection6
                    push!(c6c6, interval)
                end
                @test_throws ErrorException GenomePermutations.randomise(c6c6, f6)
                @test_nowarn GenomePermutations.randomise(c6c6, t6)
                #@test isin(GenomePermutations.randomise(GenomicFeatures.IntervalCollection([interval5]), t6))
                # can only get to 1 place if we disallow overlaps
                r1 = GenomicFeatures.IntervalCollection([
                    GenomicFeatures.Interval("chr1", 1, 80),
                    GenomicFeatures.Interval("chr1", 100, 180),
                    GenomicFeatures.Interval("chr1", 200, 280),
                    GenomicFeatures.Interval("chr1", 320, 400), 
                    GenomicFeatures.Interval("chr1", 400, 480),
                    GenomicFeatures.Interval("chr1", 500, 679)
                ])
                tr1 = GenomePermutations.StartMixture("hgTest", r1, true; max_tries = 100)
                fr1 = GenomePermutations.StartMixture("hgTest", r1, false; max_tries = 100)
                overlaps = GenomePermutations.randomise(s1, tr1)
                no_overlaps = GenomePermutations.randomise(s1, fr1)
                @test GenomePermutations.countoverlapping(no_overlaps, no_overlaps) == 6 
                @test GenomePermutations.countoverlapping(overlaps, overlaps) >= 6 
                # are chromosome by chromosome the same when they must be?
                @test GenomePermutations.randomise(collection8, f8) == collection8 == GenomePermutations.randomise(collection8, f8_by_chr)
                end
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

    @testset "[NEW!] Generic permutation test" begin 
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
            
            f6 = GenomePermutations.StartMixture("hgTest", regions6, false)
            f7 = GenomePermutations.StartMixture("hgTest", regions7, false;  by_chromosome = true)
            f8 = GenomePermutations.StartMixture("hgTest", regions8, false;  by_chromosome = true)
            f11 = GenomePermutations.StartMixture("hgTest", regions11, false;  by_chromosome = true)            
            f12 = GenomePermutations.StartMixture("hgTest", regions12, false;  by_chromosome = true)

            @testset "Single sequence" begin 
                # 1 distance to iteslf should be as significant as it gets
                test =  permtest(collection1, collection1, f6, 100) 
                @test pvalue(test.test) < 0.01
                # 2 closer than expected
                test = permtest(collection1, scollection1, f6, 100) 
                @test pvalue(test.test) < 0.01
                # 3 more distant than expected
                test = permtest(collection1, scollection3, f6, 100)
                @test pvalue(test.test) < 0.01
                # 4 sequenecs to randomise do not belong to regions 
                @test_throws Exception permtest(collection1,  collection2, f8, 100)
                # 5 Regions to randomise have an extra sequence
                @test_nowarn permtest(collection1,  collection2, f11, 100)  
            end
            @testset "Mutliple sequences" begin 
                # 1 distance to iteslf should be as significant as it gets
                test =  permtest(collection10, collection10, f7, 100) 
                @test pvalue(test.test) < 0.01
                # 2 closer than expected
                test = permtest(collection10, scollection2, f7, 100) 
                @test pvalue(test.test) < 0.01
            # 3 more distant than expected
                test = permtest(scollection5, scollection4, f7, 100)
                @test pvalue(test.test) < 0.01
                # Need a test for the case where we expect a p of 1. 
                # 4 All sequences wrong
                @test_throws Exception permtest(collection14,  collection15, f8, 100)
                # 5 One of 2 sequences are wrong
                @test_throws Exception permtest(collection14,  collection15, f11, 100)
                # 6 there is an extra chromosome - This should work fine
                @test_nowarn permtest(collection14,  collection15, f12, 100)
            end
        end
    end

    @testset "legacy functions" begin

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

    @testset "[TMP] Refactoring test" begin 
        collection1 = GenomicFeatures.IntervalCollection([
            GenomicFeatures.Interval("chr1", 1, 100),
            GenomicFeatures.Interval("chr1", 200, 300),
            GenomicFeatures.Interval("chr1", 400, 500)
            ]) # test1
        instance = GenomePermutations.StartMixture("hgTest", collection1, false)
        @test instance._distribution["chr1"] == generatedistribution(collection1)
    
    end
end
