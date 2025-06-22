def lcs_and_distance(S1, S2, d=1, r=2, e=0):
    n, m = len(S1), len(S2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]


    for i in range(n + 1): #limits
        dp[i][0] = i * d
    for j in range(m + 1): #limits
        dp[0][j] = j * d

    #dynamic table
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if S1[i - 1] == S2[j - 1]:
                cost = e  # match
            else:
                cost = r  # replace

            dp[i][j] = min(
                dp[i - 1][j] + d,        # delete
                dp[i][j - 1] + d,        # insert
                dp[i - 1][j - 1] + cost  # match or replace
            )


    distance = dp[n][m]  #distance
    u = (n + m - distance) // 2
    i, j = n, m
    lcs_str = [] #LCS recon

    while i > 0 and j > 0:
        if S1[i - 1] == S2[j - 1]:
            lcs_str.append(S1[i - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + d:
            i -= 1
        elif dp[i][j] == dp[i][j - 1] + d:
            j -= 1
        else:
            i -= 1
            j -= 1

    lcs_str.reverse()

    return distance, u, ''.join(lcs_str)


S1 = "alabama"
S2 = "alacazam"

distance, u, lcs_str = lcs_and_distance(S1, S2)

print("Edit Distance D(n, m):", distance)
print("LCS length u:", u)
print("LCS:", lcs_str)
