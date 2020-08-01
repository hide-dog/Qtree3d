using ProgressMeter
using Printf

function main()
    PARAMDAT = "PARAMDAT.json"

    # 膨張率 [%]
    expantion = 10
    
    # 構造体メモリ
    memori = Int(5.0e6)
    
    # ----------------------------------------------------------
    # 変数はないはず
    # ----------------------------------------------------------
    out_file,in_file,n_div,s_x,e_x,n_x,s_y,e_y,n_y,s_z,e_z,n_z = input_para(PARAMDAT)
    main_pre(n_x,n_y,s_x,e_x,s_y,e_y,s_z,e_z,n_z)

    nodes, nodes_num, element, element_shape = read_allgrid(in_file)

    println("nodes    : "*string(size(nodes)[1]))
    println("elements : "*string(size(element)[1]))
    

    # 空間分割
    divnum = 2^n_div

    # 辺の長さ
    lx = abs(e_x - s_x)
    ly = abs(e_y - s_y)
    lz = abs(e_z - s_z)
    lx_exp = lx *expantion/200
    ly_exp = ly *expantion/200
    lz_exp = lz *expantion/200

    # モートン単位長
    unitx = abs(e_x - s_x) / divnum *(1 + expantion/100)
    unity = abs(e_y - s_y) / divnum *(1 + expantion/100)
    unitz = abs(e_z - s_z) / divnum *(1 + expantion/100)

    # 空間分割中心算出
    dx = (e_x - s_x) / n_x
    dy = (e_y - s_y) / n_y
    dz = (e_z - s_z) / n_z

    new_morton_num = zeros(n_x,n_y,n_z,2)
    new_cell_center = zeros(n_x,n_y,n_z,3)
    for i in 1:n_x
        for j in 1:n_y
            for k in 1:n_z
                # 左上の点のモートン空間番号
                ulx = floor((s_x + dx * (i-1) -s_x + lx_exp) / unitx) # 原点に合わせて番号を付けるため（-s_x + lx_exp）
                uly = floor((s_y + dy *  j    -s_y + ly_exp) / unity)
                ulz = floor((s_z + dz *  k    -s_z + lz_exp) / unitz)

                # 右下の点のモートン空間番号
                lrx = floor((s_x + dx *  i    -s_x + lx_exp) / unitx) # 原点に合わせて番号を付けるため（-s_x + lx_exp）
                lry = floor((s_y + dy * (j-1) -s_y + ly_exp) / unity)
                lrz = floor((s_z + dz * (k-1) -s_z + lz_exp) / unitz)

                # 0<=s, 0<=t
                #new_morton_num[i,j,1],new_morton_num[i,j,2] = cal_morton(ulx,uly,lrx,lry)
                s,t = cal_morton(ulx,uly,ulz,lrx,lry,lrz,n_div)
                new_morton_num[i,j,k,1] = s
                new_morton_num[i,j,k,2] = t

                new_cell_center[i,j,k,1] = s_x +dx/2+dx*(i-1)
                new_cell_center[i,j,k,2] = s_y +dy/2+dy*(j-1)
                new_cell_center[i,j,k,3] = s_z +dz/2+dz*(k-1)
            end
        end
    end
    println(" fin kuukan bunkatu\n")

    # 登録管理構造体
    s = Int64(0)
    for i in 1:(n_div+1)
        s = s + 8^(i-1)
    end
    morton_tree = zeros(Int64,s,memori)
    morton_tree_num = ones(Int64,s)

    # セルを囲む直方体
    # これでモートン空間番号を算出
    cellnum = size(element)[1]
    cell_center = zeros(cellnum,3)
    
    for i in 1:cellnum
        loop = element_shape[i]
        x = zeros(loop)
        y = zeros(loop)
        z = zeros(loop)

        for j in 1:loop
            x[j] = nodes[element[i,j],1]
            y[j] = nodes[element[i,j],2]
            z[j] = nodes[element[i,j],3]

            cell_center[i,1] += x[j] / loop
            cell_center[i,2] += y[j] / loop
            cell_center[i,3] += z[j] / loop
        end
        
        maxx = maximum(x)
        minx = minimum(x)
        maxy = maximum(y)
        miny = minimum(y)
        maxz = maximum(z)
        minz = minimum(z)

        # 左上の点のモートン空間番号
        ulx = floor((minx-s_x + lx_exp) / unitx) # 原点に合わせて番号を付けるため（-s_x + lx_exp）
        uly = floor((maxy-s_y + ly_exp) / unity)
        ulz = floor((maxz-s_z + lz_exp) / unitz)

        # 右下の点のモートン空間番号
        lrx = floor((maxx-s_x + lx_exp) / unitx) # 原点に合わせて番号を付けるため（-s_x + lx_exp）
        lry = floor((miny-s_y + ly_exp) / unity)
        lrz = floor((minz-s_z + lz_exp) / unitz)

        # =-------------異常値をはじく事-----------
        if (ulx < 0 || lrx < 0) || (uly < 0 || lry < 0) || (ulz < 0 || lrz < 0)
            continue
        elseif (ulx > divnum || lrx > divnum) || (uly > divnum || lry > divnum) || (ulz > divnum || lrz > divnum)
            continue
        else
            #morton_num[i,1],morton_num[i,2] = cal_morton(ulx,uly,lrx,lry)
            # 0<=s, 0<=t
            s,t = cal_morton(ulx,uly,ulz,lrx,lry,lrz,n_div)

            num = Int64(t + (8^s-1)/7   +1)     # 線形空間番号

            s = morton_tree_num[num]
            morton_tree[num,s] = i
            morton_tree_num[num] += 1
        end
    end
    println(" fin kuukan bunkatu2 \n")
    
        
    # 線形探索用
    s,search0,search1,search2,search3,search4 = liner_morton(n_div)

    # 分割したセル中心と近い点を四点計算
    point = 4
    close_point_num = zeros(n_x,n_y,n_z,point)
    distance = zeros(n_x,n_y,n_z,point) 
    distance[:,:,:,:] .= 100
    
    prog = Progress(n_x,1)
    for i in 1:n_x
        next!(prog)
        for j in 1:n_y
            for k in 1:n_z
                if new_morton_num[i,j,k,1] ==0      # 所属空間
                    temp = copy(search0)
                    loop = Int(s[1] * 8^0)
                elseif new_morton_num[i,j,k,1] ==1      # 所属空間
                    temp = copy(search1)
                    loop = Int(s[2] * 8^1)
                elseif new_morton_num[i,j,k,1] ==2      # 所属空間
                    temp = copy(search2)
                    loop = Int(s[3] * 8^2)
                elseif new_morton_num[i,j,k,1] ==3      # 所属空間
                    temp = copy(search3)
                    loop = Int(s[4] * 8^3)
                elseif new_morton_num[i,j,k,1] ==4      # 所属空間
                    temp = copy(search4)
                    loop = Int(s[5] * 8^4)
                else 
                    "\n check new_morton programm \n"
                    throw(UndefVarError(:x))
                end

                # 空間事に探査
                for pointi in 1:loop
                    for l in 1:memori
                        
                        if morton_tree[temp[pointi],l] == 0
                            break
                        end
                        bangou = morton_tree[temp[pointi],l]

                        x1 = new_cell_center[i,j,k,1]
                        y1 = new_cell_center[i,j,k,2]
                        z1 = new_cell_center[i,j,k,3]
                        x2 = cell_center[bangou,1]
                        y2 = cell_center[bangou,2]
                        z2 = cell_center[bangou,3]
                        temp_distance = ((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 )^0.5

                        for ii in 1:point
                            if temp_distance < distance[i,j,k,ii]
                                for m in 1:(point+1-ii)-1
                                    distance[i,j,k,point+1-m] = distance[i,j,k,point-m]
                                    close_point_num[i,j,k,point+1-m] = close_point_num[i,j,k,point-m]
                                end
                                distance[i,j,k,ii] = temp_distance
                                close_point_num[i,j,k,ii] = morton_tree[temp[pointi],l]
                                break
                            end
                        end
                    end
                end
            end
        end
    end

    println(" fin kuukan bunkatu3 \n")
    
    #=
    println("\n ---------------------------------- \n")
    println(" start interpolation ")
    println("\n ---------------------------------- \n")
    
    # 補間 loop
    Qcell = zeros(n_x,n_y,5)
    for i in 1:n_x
        for j in 1:n_y
            # cell中心ベクトル p
            p  = [new_cell_center[i,j,1], new_cell_center[i,j,2]]
            
            # 四点ベクトル p
            x0 = cell_center[Int64(close_point_num[i,j,1]),1]
            y0 = cell_center[Int64(close_point_num[i,j,1]),2]
            x1 = cell_center[Int64(close_point_num[i,j,2]),1]
            y1 = cell_center[Int64(close_point_num[i,j,2]),2]
            x2 = cell_center[Int64(close_point_num[i,j,3]),1]
            y2 = cell_center[Int64(close_point_num[i,j,3]),2]

            p0 = [x0, y0]
            p1 = [x1, y1]
            p2 = [x2, y2]
            
            # 局所座標ベクトル e
            eu = [p1[1] - p0[1], p1[2] - p0[2]]
            ev = [p2[1] - p0[1], p2[2] - p0[2]]
            inv_A = inverse_matrix(eu,ev)

            # 補間距離ベクトル u
            up = zeros(2)
            for l in 1:2
                up[l] = inv_A[l,1]*(p[1]-p0[1]) + inv_A[l,2]*(p[2]-p0[2])
            end
            for k in 1:5
                Qcell[i,j,k] = readQ[Int64(close_point_num[i,j,1]),k]
                            + up[1]*readQ[Int64(close_point_num[i,j,2]),k]
                            + up[2]*readQ[Int64(close_point_num[i,j,3]),k]
                #Qcell[i,j] = Qcell[i,j]*avoga
            end
        end
    end
    =#

    println("\n ---------------------------------- \n")
    println(" start atari ")
    println("\n ---------------------------------- \n")

    #println(nodes)
    #println(close_point_num)

    # 分割した空間がセルと当たるか否かの判定
    atari = zeros(Int64,n_x,n_y,n_z)
    
    for i in 1:n_x
        for j in 1:n_y
            for k in 1:n_z
                for l in 1:point
                    # ---------- 過大評価 ---------------
                    cnum = Int64(close_point_num[i,j,k,l])

                    loop = element_shape[cnum]
                    x = zeros(loop)
                    y = zeros(loop)
                    z = zeros(loop)
                    
                    for m in 1:loop
                        x[m] = nodes[element[cnum,m],1]
                        y[m] = nodes[element[cnum,m],2]
                        z[m] = nodes[element[cnum,m],3]
                    end

                    maxx = maximum(x)
                    minx = minimum(x)
                    maxy = maximum(y)
                    miny = minimum(y)
                    maxz = maximum(z)
                    minz = minimum(z)
                    
                    x = new_cell_center[i,j,k,1]
                    y = new_cell_center[i,j,k,2]
                    z = new_cell_center[i,j,k,3]
                    

                    if minx < x && x < maxx
                        if miny < y && y < maxy
                            if minz < z && z < maxz
                                atari[i,j,k] += 1
                            end
                        end
                    end
                end
            end
        end
    end

    # 修正当たり判定
    #=
    atari_temp = zeros(n_x,n_y,n_z)
    for i in 2:n_x-1
        for j in 2:n_y-1
            for k in 2:n_z-1
                if atari[i,j,k] == 0 && atari[i+1,j,k] ==1 && atari[i,j+1,k] ==1 && atari[i-1,j,k] ==1 && atari[i,j-1,k] ==1 && atari[i,j,k+1] ==1 && atari[i,j,k-1] ==1
                    atari_temp[i,j,k] = 1
                else atari[i,j,k] == 1 && atari[i+1,j,k] ==0 && atari[i,j+1,k] ==0 && atari[i-1,j,k] ==0 && atari[i,j-1,k] ==0 && atari[i,j,k+1] ==0 && atari[i,j,k-1] ==0
                    atari_temp[i,j,k] = 0
                end
            end
        end
    end
    
    for i in 2:n_x-1
        for j in 2:n_y-1
            for k in 2:n_z-1
                atari[i,j,k] = atari[i,j,k] * atari_temp[i,j,k]
            end
        end
    end
    =#



    # output
    fff = "grid_morton/atari"
    open(fff,"w") do f
        write(f,"result:atari\n")
        for i in 1:n_x
            for j in 1:n_y
                for k in 1:n_z
                    a1 = string(atari[i,j,k])
                    write(f, a1*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
    
    # ダミー
    Qcell = ones(n_x,n_y,n_z,5)
    
    fff = "result_morton/"* out_file
    open(fff,"w") do f
        write(f,"result:rho[kg/m^3], u[m/s], v[m/s], p[Pa], T[K]\n")
        for i in 1:n_x
            for j in 1:n_y
                for k in 1:n_z
    
                    a1 = @sprintf("%8.8f", Qcell[i,j,k,1])
                    a2 = @sprintf("%8.8f", Qcell[i,j,k,2])
                    a3 = @sprintf("%8.8f", Qcell[i,j,k,3])
                    a4 = @sprintf("%8.8f", Qcell[i,j,k,4])
                    a5 = @sprintf("%8.8f", Qcell[i,j,k,5])

                    #=
                    a1 = string(Int64(close_point_num[i,j,1]))
                    a2 = string(Int64(close_point_num[i,j,2]))
                    a3 = string(Int64(close_point_num[i,j,3]))
                    =#
        
                    write(f, a1*" "*a2*" "*a3*" "*a4* " "*a5*"\n")
                end
            end
        end
    end
    println("\nwrite "*fff)
    
        
end



# ----------------
main()
# ----------------
