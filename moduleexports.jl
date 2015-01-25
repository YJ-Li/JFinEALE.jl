

module A

module a
function test2(n::String)
   println("$n in module A")
end
export test2
end    # end module C
using A.a
export test2

# end

# module B

module c
import A.test2
function test2(n::Float64)
   println("$n in module c")
end
export test2
end    # end module C
using A.c
export test2

module b
import A.test2
function test2(n::Int64)
   println("$n in module b")
end
export test2
end    # end module C
using A.b
export test2

end
