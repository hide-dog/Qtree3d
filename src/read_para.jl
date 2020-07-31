import JSON

function make_json(PARAMDAT)
    fff=[]
    open(PARAMDAT, "r") do f
        fff=read(f,String)
    end

    fff=replace(fff,r"#(.+)\n" => "\n")
    
    read_PARAMDAT="read_"*PARAMDAT
    open(read_PARAMDAT, "w") do f
        write(f,fff)
    end
    return read_PARAMDAT
end 

function read_json(read_PARAMDAT)
    dict=1
    open(read_PARAMDAT, "r") do f
        dicttxt = read(f,String)  # file information to string
        dict=JSON.parse(dicttxt)  # parse and transform data
    end
    return dict
end

function read_para(dict)
    # -- file --
    out_file = dict["out_file_front"]*dict["out_ext"]
    in_file = dict["in_file_data"]
    
    # -- 分割数　--
    n_div = Int(parse(Float64,dict["n_div"]))
    
    # -- 空間分割 --
    s_x = parse(Float64,dict["s_x"])
    e_x = parse(Float64,dict["e_x"])
    n_x = Int(parse(Float64,dict["n_x"]))

    s_y = parse(Float64,dict["s_y"])
    e_y = parse(Float64,dict["e_y"])
    n_y = Int(parse(Float64,dict["n_y"]))

    s_z = parse(Float64,dict["s_z"])
    e_z = parse(Float64,dict["e_z"])
    n_z = Int(parse(Float64,dict["n_z"]))
    

    return out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y,s_z,e_z,n_z
end

function input_para(PARAMDAT)
   read_PARAMDAT=make_json(PARAMDAT)
   dict=read_json(read_PARAMDAT)
   out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y,s_z,e_z,n_z = read_para(dict)
   println("fin read para")
   return out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y,s_z,e_z,n_z
end