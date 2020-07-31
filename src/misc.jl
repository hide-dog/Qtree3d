function ten_to_two(n)
    b = zeros(Int8,15)
    temp_b = zeros(Int8,15)

    ite=1
    while n > 0
        b[ite] = n % 2
        n = Int((n-b[ite])/2)
        ite = ite +1
    end
    
    #=
    for i in 1:6
        temp_b[i] = b[6+1-i]
    end
    =#
    
    return b
end

function two_to_ten(b)
    n = Int8(0)

    for i in 1:length(b)
        n = n + 2^(i-1)*b[i]
    end
    return n
end

# ----------------------
# -- Inverse matrix   --
# ----------------------
function inverse_matrix(eu,ev)
# ----------------------
# -- A = ( eu[1]  ev[1])  --
# --     ( eu[2]  ev[2])  --
# ----------------------
    invA = zeros(2,2)
    detA = eu[1]*ev[2] - eu[2]*ev[1]

    invA[1,1] = (ev[2])/detA
    invA[1,2] = (-ev[1])/detA
    invA[2,1] = (-eu[2])/detA
    invA[2,2] = (eu[1])/detA
    return invA
end
        

function cal_morton(ulx,uly,ulz,lrx,lry,lrz,n_div) # 0<=s,tでreturn
    # 左上奥と右下前の座標を整数で
    morton_belong = 0
    morton_num = 0

    upl_morton = zeros(Int8,15)
    lowr_morton = zeros(Int8,15)
    temp_morton = zeros(Int8,15)


    # 左上の点のモートン空間番号
    x_two = ten_to_two(ulx)
    y_two = ten_to_two(uly)
    z_two = ten_to_two(ulz)
    
    for k in 1:5
        upl_morton[k*3-2] = x_two[k]
        upl_morton[k*3-1] = y_two[k]
        upl_morton[k*3]   = z_two[k]
    end

    # 右下の点のモートン空間番号
    x_two = ten_to_two(lrx)
    y_two = ten_to_two(lry)
    z_two = ten_to_two(lrz)
    
    for k in 1:5
        lowr_morton[k*3-2] = x_two[k]
        lowr_morton[k*3-1] = y_two[k]
        lowr_morton[k*3]   = z_two[k]
    end
    
    # 排他的論理和 xor
    for i in 1:15
        if upl_morton[i]+lowr_morton[i] == 2
            temp_morton[i] = 0
        else
            temp_morton[i] = upl_morton[i]+lowr_morton[i]
        end
    end
    
    # 所属空間のチェック
    
    ite = Int64(0)
    for i in 1:n_div
        s30 = 3*(n_div-(i-1))
        s31 = 3*(n_div-(i-1)) - 1
        s32 = 3*(n_div-(i-1)) - 2
        
        if temp_morton[s30] == 1 || temp_morton[s31] == 1 || temp_morton[s32] == 1
            morton_belong = ite                             # ルート空間
            break
        end
        #=
        if temp_morton[s_even] == 1 || temp_morton[s_odd] == 1
            morton_belong = 0                             # ルート空間
        elseif temp_morton[4] == 1 || temp_morton[3] == 1
            morton_belong = 1                             # 親空間
        elseif temp_morton[2] == 1 || temp_morton[1] == 1
            morton_belong = 2                             # 子空間
        else
            morton_belong = 3                             # 孫空間
        end
        =#
        ite += 1
    end
    
    # 所属空間番号
    temp_morton = zeros(Int8,15)

    s30 = 3*(n_div)
    s31 = 3*(n_div) - 1
    s32 = 3*(n_div) - 2
    if morton_belong == 0                             # ルート空間
        morton_num = 0
    elseif morton_belong == 1                         # 親空間
        # 4右シフト
        temp_morton[1] = lowr_morton[s32]
        temp_morton[2] = lowr_morton[s31]
        temp_morton[3] = lowr_morton[s30]
        
        morton_num = two_to_ten(temp_morton)

    elseif morton_belong == 2                         # 子空間        
        # 2右シフト
        temp_morton[1] = lowr_morton[s32-3]
        temp_morton[2] = lowr_morton[s31-3]
        temp_morton[3] = lowr_morton[s30-3]
        temp_morton[4] = lowr_morton[s32]
        temp_morton[5] = lowr_morton[s31]
        temp_morton[6] = lowr_morton[s30]
        
        morton_num = two_to_ten(temp_morton)
    elseif morton_belong == 3                         # 孫空間
        
        # 0右シフト
        temp_morton[1] = lowr_morton[s32-6]
        temp_morton[2] = lowr_morton[s31-6]
        temp_morton[3] = lowr_morton[s30-6]
        temp_morton[4] = lowr_morton[s32-3]
        temp_morton[5] = lowr_morton[s31-3]
        temp_morton[6] = lowr_morton[s30-3]
        temp_morton[7] = lowr_morton[s32]
        temp_morton[8] = lowr_morton[s31]
        temp_morton[9] = lowr_morton[s30]
        
        morton_num = two_to_ten(temp_morton)
        
    elseif morton_belong == 4                         # 孫空間
        temp_morton = copy(lowr_morton)
        morton_num = two_to_ten(temp_morton)
    end

    return morton_belong,morton_num
end


function liner_morton(n_div)
    # 両方の探査を行うため，煩雑
    s = zeros(Int64,5)

    for j in 1:5
        for i in 1:(n_div+1)        
            s[j] += floor(8.0^(i-j))
        end
        s[j] += (j-1)
    end

    search0 = zeros(Int64,8^0,s[1])
    search1 = zeros(Int64,8^1,s[2])
    search2 = zeros(Int64,8^2,s[3])
    search3 = zeros(Int64,8^3,s[4])
    search4 = zeros(Int64,8^4,s[5])
    
    # ルート空間：s0
    for i in 1:1
        for j in 1:s[1]
            search0[i,j] = j
        end
    end

    # 親空間：s1
    for i in 1:8
        search1[i,1] = 1　　　　　　　　　　　　　　　　# ルート空間
        search1[i,2] = i+1                           # 親空間
        for j in 1:8
            search1[i,j+2] = 10 + 8*(i-1) + (j-1)     # 子空間
        end
        for j in 1:64
            search1[i,j+10] = 73 + 8^2*(i-1) + (j-1)   # 孫空間
        end

        if n_div ==4
            for j in 1:512
                search1[i,j+74] = 1+8+64+512+1 + 8^3*(i-1) + (j-1)   # ひ孫空間
            end
        end
    end
    
    # 子空間：s2
    for i in 1:64
        search2[i,1] = 1                             # ルート空間
        search2[i,2] = div(i-1,8) + 2                # 親空間
        search2[i,3] = i+9                           # 子空間
        for j in 1:8
            search2[i,j+3] = 10 + 8*(i-1) + (j-1)    # 孫空間
        end
        if n_div ==4
            for j in 1:64
                search2[i,j+11] = 73 + 8^2*(i-1) + (j-1)    # 孫空間
            end
        end
    end
    
    # 孫空間
    for i in 1:512
        search3[i,1] = 1                             # ルート空間
        search3[i,2] = div(i-1,64) + 2               # 親空間
        search3[i,3] = div(i-1,8) + 10                # 子空間
        search3[i,4] = 73 + (i-1)                    # 孫空間
        if n_div >=4
            for j in 1:8
                search3[i,j+4] = 10 + 8*(i-1) + (j-1)    # 孫空間
            end
        end
    end
    
    # ひ孫空間
    if n_div >=4
        for i in 1:4096
            search4[i,1] = 1                             # ルート空間
            search4[i,2] = div(i-1,512) + 2               # 親空間
            search4[i,3] = div(i-1,64) + 10                # 子空間
            search4[i,4] = div(i-1,8) + 73               # 孫空間
            search4[i,5] = 585 + (i-1)                    # 孫空間
        end
    end

    return s,search0,search1,search2,search3,search4
end
