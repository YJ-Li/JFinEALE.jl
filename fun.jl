function test1(f, n)
    for i=1:n
        f(i)
    end
end
function test2(n)
    for i=1:n
        g(i)
    end
end
function test3(f, n)
    for i=1:n
        invoke(f, (Int,), i)
    end
end
g(i) = i


function test4{T}(o::T, n)
    for i=1:n
        get(o,i)
    end
end
type MyFun1
end
function get(o::MyFun1,i::Int)
    return g(i)
end
function test5{T}(o::T, n)
    r=0
    for i=1:n
        update!(o,r,i)
    end
end
type MyFun2
end
function update!(o::MyFun2,result::Int,i::Int)
    result=g(i)
    return 
end
n=10000000;
@time test1(g, n)
@time test2(n)
@time test3(g, n)
@time test4(MyFun1(), n)
@time test5(MyFun2(), n)
     
