using Nonconvex: is_array_of_points

@testset "Array utils" begin
    @test is_array_of_points([[1], [2], [3]]) == true
    @test is_array_of_points([1, 2, 3, 4, 5, 6]) == true
    @test is_array_of_points([[1,2], [3,4], [4,5]]) == true
    @test is_array_of_points([[[1, 2], [3, 4], [5, 6]], [[1, 2], [3, 4], [5, 6]], 
                    [[1, 2], [3, 4], [5, 6]], [[1, 2], [3, 4], [5, 6]]]) == false
    @test is_array_of_points([[[1, 2], [3, 4], [5, 6]], [[1, 2], [3, 4], [5, 6]], 
                    [[1, 2], [3, 4], [5, 6]], [[1, 2], [3, 4], [5, 6, 7]]]) == false
    @test is_array_of_points([[1, 2, 3], [4, 5], [6, 7, 8]]) == false
    @test is_array_of_points(["a", 2, 3]) == false
end
