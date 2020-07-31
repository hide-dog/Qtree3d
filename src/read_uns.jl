# ----------------------
# -- read             --
# ----------------------

function read_uns(skipnum,uns)
    # skipnum = 9
    ite = Int64(0)

    fff=[]
    open(uns, "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum

    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
        fff[i]=replace(fff[i],"\r" => "")
    end

    ite += skipnum

    # boundary
    ite += 1
    temp = split(fff[ite]," ")
    bdnum = parse(Int64,temp[1])
    for i in 1:bdnum
        # 境界条件処理
        ite += 1
    end
    
    # skip
    ite += 4

    # Nodes
    ite += 1
    temp = split(fff[ite]," ")
    nodes_num = parse(Int64,temp[1])
    nodes = zeros(nodes_num,3)

    for i in 1:nodes_num
        ite +=1
        temp = split(fff[ite]," ")

        k = 1
        for j in 1:length(temp)
            if temp[j] != ""
                nodes[i,k] = parse(Float64,temp[j])
                k += 1
            end
        end
    end

    println("read nodes ")
    
    #=
    k = 1
    for j in 1:length(temp)
        if temp[j] != ""
            readQ[i,k] = parse(Float64,temp[j])
            k += 1
        end
    end
    =#

    # boundary face
    ite += 2
    temp = split(fff[ite]," ")
    bdface_num = parse(Int64,temp[1])
    for i in 1:bdface_num
        # 境界条件処理
        ite += 1
    end

    # elements
    ite += 1
    temp = split(fff[ite]," ")
    
    k=0
    ite_temp = ite
    while k == 0
        if occursin("Variables",fff[ite_temp]) == true
            k=1
        end
        
        ite_temp += 1
    end

    println("read elements ")

    ele_num = ite_temp - ite -2
    element = zeros(Int64,ele_num,8)
    element_shape = zeros(Int64,ele_num)

    for i in 1:ele_num
        ite +=1
        temp = split(fff[ite]," ")

        nnn = parse(Int64,temp[1])
        if nnn == 1
            element_shape[i] = 4
        elseif nnn == 2
            element_shape[i] = 8
        elseif nnn == 3
            element_shape[i] = 6
        elseif nnn == 4
            element_shape[i] = 5
        else
            " check element in .uns !! \n"
            throw(UndefVarError(:x))
        end
        
        for j in 1:element_shape[i]
            element[i,j] = parse(Int64,temp[2+j])
        end
    end
    
    return nodes, nodes_num, element, element_shape
end


function read_result(skipnum)
    fff=[]
    open("test333.dat", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_point=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    readQ = zeros(num_point,5)
    for i in 1:num_point
        temp=split(fff[i+skipnum]," ")
        
        k = 1
        for j in 1:length(temp)
            if temp[j] != ""
                readQ[i,k] = parse(Float64,temp[j])
                k += 1
            end
        end

    end
    return readQ
end

function read_allgrid(uns)
    skip=9
    nodes, nodes_num, element, element_shape = read_uns(skip,uns)

    #readQ     = read_result(skip)
    println("fin read grid")
    return nodes, nodes_num, element, element_shape
end