$files = Get-ChildItem "C:\Users\Salazar\Documents\tesis\tesis\instances\250_40_30"
foreach ($f in $files){
    $content = $f.FullName
    Write-Output $content
    Invoke-Expression -Command "julia.exe --sysimage JuliaSysImage.dll src\heuristics\cons_local.jl $content "

}