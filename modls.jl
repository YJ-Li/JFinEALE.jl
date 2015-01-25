

module A

module B
function test2(n)
   println("$n in module B")
end
export test2
end   # end module B

module C
function test3(n)
   println("$n in module C")
end
export test3
end    # end module C

using .B
export test2
using A.C
export test3

end     # end module A

module AZ

module B
a=1
export a
end

using AZ.B
module C
using AZ.B
b=a+1
end

end
