using Printf

function main_pre(n_x,n_y,s_x,e_x,s_y,e_y,s_z,e_z,n_z)
    xy_or_r = 1          # 1:x,y 2:r,t
    xnum = n_x           #       2:rnum
    ynum = n_y           #       2:tnum
    znum = n_z           #       2:tnum
    lenx = abs(e_x-s_x)           #       2:inr
    leny = abs(e_y-s_y)           #       2:outr
    lenz = abs(e_z-s_z)           #       2:outr
    st_ang  = -1/2*pi       # 2:angle
    en_ang  = -3/2*pi          # 2:angle

    outdir = "grid_morton"
    make_dir(outdir)
    result = "result_morton"
    make_dir(result)
    result = "post_result_morton"
    make_dir(result)
    nodes, xnum_max, ynum_max, znum_max = mk_gird(xnum,ynum,znum,lenx,leny,lenz,xy_or_r,s_x,s_y,s_z,st_ang,en_ang,outdir)
end

function mk_gird(xnum,ynum,znum,lenx,leny,lenz,xy_or_r,s_x,s_y,s_z,st_ang,en_ang,outdir)
    """
    nodes[i,j,k]
    i   : x,r方向の番号
    j   : y,theta方向の番号
    k=1 : 点のx座標
    k=2 : 点のy座標

    nodes[1,:]      : x方向境界
    nodes[:,1]      : y方向境界
    nodes[xnum+2,:] : x方向境界
    nodes[:,ynum+2] : y方向境界
    """
    xnum_max = xnum+1+2
    ynum_max = ynum+1+2
    znum_max = znum+1+2
    nodes = zeros(xnum_max, ynum_max, znum_max,3)
    if xy_or_r == 1
        for i in 2:xnum_max-1
            for j in 2:ynum_max-1
                for k in 2:znum_max-1
                    x = lenx/(xnum)*(i-2) + s_x
                    y = leny/(ynum)*(j-2) + s_y
                    z = lenz/(znum)*(k-2) + s_z
                    nodes[i,j,k,1] = x
                    nodes[i,j,k,2] = y
                    nodes[i,j,k,3] = z
                end
            end 
        end
    
    end

    #=
    仮想セルの作成：現在あるセル境界線を延長して作成
    そのまま延長した際にセルが交差するのを防ぐため、仮想セルの大きさを100^-1倍にする.
    ベクトルで考えれば下記のような計算になる．たぶん
    また、[1,1]等の角にはダミー値が入っている.
    =#
    #=
    for j in 1:ynum_max
        for k in 1:2
            nodes[1,j,k]        =  1.01*nodes[2,j,k] - 0.01*nodes[3,j,k]
            nodes[xnum_max,j,k] =  1.01*nodes[xnum_max-1,j,k] - 0.01*nodes[xnum_max-2,j,k]
        end
    end
    for i in 1:xnum_max
        for k in 1:2
            nodes[i,1,k]        =  1.01*nodes[i,2,k] - 0.01*nodes[i,3,k]
            nodes[i,ynum_max,k] =  1.01*nodes[i,ynum_max-1,k] - 0.01*nodes[i,ynum_max-2,k]
        end
    end
    =#
    
    #=
    fff=outdir*"/nodes"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum , x, y\n")
        for i in 1:xnum_max
            for j in 1:ynum_max
                x = @sprintf("%8.8e", nodes[i,j,1])
                y = @sprintf("%8.8e", nodes[i,j,2])
                write(f,string(i)*" "*string(j)*" "*x*" "*y*"\n")
            end
        end
    end
    println("write "*fff)
    =#

    # nodes_forvtk
    nodes_num = zeros(Int,xnum_max, ynum_max, znum_max)
    fff=outdir*"/nodes_forvtk"
    open(fff,"w") do f
        write(f,"nodes: xnum, ynum , x, y\n")
        a=1
        for i in 2:xnum_max-1
            for j in 2:ynum_max-1
                for k in 2:znum_max-1
                    x = @sprintf("%8.8e", nodes[i,j,k,1])
                    y = @sprintf("%8.8e", nodes[i,j,k,2])
                    z = @sprintf("%8.8e", nodes[i,j,k,3])
                    write(f,string(a)*" "*x*" "*y*" "*z*"\n")
                    nodes_num[i,j,k] = a
                    a = a+1
                end
            end
        end
    end
    println("write "*fff)

    fff=outdir*"/nodesnum"
    open(fff,"w") do f
        write(f,"nodesnum: xnum_max, ynum_max\n")
        write(f,string(xnum_max)*" "*string(ynum_max)*" "*string(znum_max)*"\n")
    end

    ### element ###
    fff=outdir*"/element_forvtk"
    open(fff,"w") do f
        write(f,"elements:cell_xnum, l1,l2,l3,l4,u1,u2,u3,u4 \n")
        a=1
        for i in 2:xnum_max-2
            for j in 2:ynum_max-2
                for k in 2:znum_max-2
                    write(f,string(a))
                    for l in 1:2
                        d1 = @sprintf("%1.0f", nodes_num[i  ,j  ,k+(l-1)])
                        d2 = @sprintf("%1.0f", nodes_num[i  ,j+1,k+(l-1)])
                        d3 = @sprintf("%1.0f", nodes_num[i+1,j+1,k+(l-1)])
                        d4 = @sprintf("%1.0f", nodes_num[i+1,j  ,k+(l-1)])
                        write(f," "*d1*" "*d2*" "*d3*" "*d4)
                    end
                    write(f,"\n")
                    a = a+1
                end
            end
        end
    end
    println("write "*fff)

    return  nodes,xnum_max,ynum_max,znum_max
end

function make_dir(outdir)
    k = 0
    try rm(outdir,recursive=true)
    catch
        mkdir(outdir)
        k = 1
    end

    if k == 0
        mkdir(outdir)
    end
end
