using Printf

function main()
    output_dir_name="post_result_morton"
    input_result_dir="result_morton"

    #vtk_頭の文字
    front="# vtk DataFile Version 2.0"*"\n"
    back="\n"*"ASCII"*"\n"*"DATASET UNSTRUCTURED_GRID"*"\n"
    write_file(output_dir_name,input_result_dir,front,back)
end
    
function write_file(outdir,inresult_dir,front,back)

    nodes_xmax, nodes_ymax, nodes_zmax = read_nodenum(1)
    nodes = read_nodes(1)

    elements = read_elements(1)

    atari = read_atari(1)
    atarinum = 0
    for i in 1:length(atari)
        if atari[i] > 0
            atarinum += 1
        end
    end

    inf = readdir(inresult_dir)
    cellnum = (nodes_xmax-3)*(nodes_ymax-3)*(nodes_zmax-3)
    rho = zeros(cellnum)
    u = zeros(cellnum)
    v = zeros(cellnum)
    p = zeros(cellnum)
    T = zeros(cellnum)
    
    for i in 1:length(inf)
        if occursin(".dat", inf[i]) == true
            dname = replace(inf[i],".dat" => "")
            out_file = replace(inf[i],".dat" => ".vtk")
            print("start writing "*out_file*"\n")

            fff=[]
            open(inresult_dir*"/"*inf[i], "r") do f
                fff=read(f,String)
            end 
            fff=split(fff,"\n",keepempty=false)

            for j in 2:length(fff)
                fff[j]=replace(fff[j]," \r" => "")
                temp = split(fff[j]," ")
                rho[j-1] = parse(Float64,temp[1])
                u[j-1] = parse(Float64,temp[2])
                v[j-1] = parse(Float64,temp[3])
                p[j-1] = parse(Float64,temp[4])
                T[j-1] = parse(Float64,temp[5])
            end
            
            write_points(nodes,out_file,outdir,dname,front,back)
            write_cells(elements,atari,out_file,outdir,atarinum)
            write_result(rho,atari,out_file,outdir,1,atarinum)
            write_result(u,atari,out_file,outdir,2,atarinum)
            write_result(v,atari,out_file,outdir,3,atarinum)
            write_result(p,atari,out_file,outdir,4,atarinum)
            write_result(T,atari,out_file,outdir,5,atarinum)
            print("fin writing "*out_file*"\n")
        end
    end
end 

function read_nodenum(skipnum)
    # xmax : 仮想セルも含めたnodeのxの数
    # ymax : 仮想セルも含めたnodeのyの数    
    fff=[]
    open("grid_morton/nodesnum", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum

    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    temp = split(fff[2]," ")
    xmax = parse(Int64,temp[1]) 
    ymax = parse(Int64,temp[2]) 
    zmax = parse(Int64,temp[3])
    
    return xmax, ymax, zmax
end

function read_nodes(skipnum)
    
    fff=[]
    open("grid_morton/nodes_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_nodes=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    nodes=zeros(num_nodes,3)
    for i in 1:num_nodes
        temp=split(fff[i+skipnum]," ")

        # x = parse(Float64,temp[1])
        # y = parse(Float64,temp[2])
        # z = parse(Float64,temp[3])
        x = parse(Float64,temp[2])
        y = parse(Float64,temp[3])
        z = parse(Float64,temp[4])
        nodes[i,1] = x
        nodes[i,2] = y
        nodes[i,3] = z
    end
    return nodes
end 

function read_elements(skipnum)
    dim = 3
    a = 2^dim

    fff=[]
    open("grid_morton/element_forvtk", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num_elements=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    elements = zeros(num_elements,a)
    for i in 1:num_elements
        temp=split(fff[i+skipnum]," ")
        for j in 1:a
            elements[i,j] = parse(Int64,temp[j+1])-1
        end
    end
    return elements
end 

function read_atari(skipnum)
    fff=[]
    open("grid_morton/atari", "r") do f
        fff=read(f,String)
    end 
    fff=split(fff,"\n",keepempty=false)
    num=length(fff)-skipnum
    
    for i in 1+skipnum:length(fff)
        fff[i]=replace(fff[i]," \r" => "")
    end

    atari = zeros(num)
    for i in 1:num
        temp=split(fff[i+skipnum]," ")
        atari[i,1] = parse(Int64,temp[1])
    end
    return atari
end 


function  write_points(nodes,out_file,outdir,dname,front,back)
    a = size(nodes)[1]
    a_st = @sprintf("%1.0f", a)

    fff = outdir*"/"*out_file
    open(fff,"w") do f
        write(f,front)
        write(f,dname)
        write(f,back)
        write(f,"POINTS "*a_st*" float\n")
        for i in 1:a
            x = @sprintf("%7.7f", nodes[i,1])
            y = @sprintf("%7.7f", nodes[i,2])
            z = @sprintf("%7.7f", nodes[i,3])
            write(f, x*" "*y*" "*z*"\n")
        end
    end
    println("fin writing points")
end

function write_cells(elements,atari,out_file,outdir,atarinum)
    a = size(elements)[1]
    b = a*9                                     # dim !!!
    a_st = @sprintf("%1.0f", atarinum)
    b_st = @sprintf("%1.0f", atarinum*9)         # dim !!!

    dim = 3
    dimnum = 2^dim
    dimnum1 = 12

    fff = outdir*"/"*out_file
    open(fff,"a") do f
        write(f, "CELLS "*a_st*" "*b_st*"\n")
        for i in 1:a
            if atari[i] > 0
                write(f, string(dimnum))
                    
                for j in 1:dimnum
                    d = @sprintf("%1.0f", elements[i,j])
                    write(f," "*d)
                end
                write(f,"\n")
            end
        end
        write(f, "CELL_TYPES "*a_st*"\n")
        for i in 1:a
            if atari[i] > 0
                write(f, string(dimnum1)*"\n")          #四角のみ
            end
        end
    end     
    println("fin writing cells")
end

function write_result(val,atari,out_file,outdir,k,atarinum)
    dtype = ["rho","u","v","p","T"]

    a = length(val)
    a_st = @sprintf("%1.0f", atarinum)
    
    temp1 = "CELL_DATA "*a_st*"\n"
    temp2 = "SCALARS "*dtype[k]*" float\nLOOKUP_TABLE default\n"

    fff = outdir*"/"*out_file
    open(fff,"a") do f
        if k == 1
            write(f, temp1)
        end
        write(f, temp2)
        for i in 1:a
            if atari[i] > 0
                data = @sprintf("%7.7f", val[i])
                write(f, data)
                write(f, "\n")
            end
        end
    end 
end
# main
main()
