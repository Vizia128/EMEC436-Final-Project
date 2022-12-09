using Flux, JLD2, Random, LinearAlgebra, Statistics, Plots
using Flux.Optimise

device = cpu

function load_data()
    @load "training_data/airfoil_training_data.jld2" x_sdf y_puv

    split = Int(floor(0.8*size(y_puv, 4)))
    train_data = Flux.DataLoader(((x_sdf[1][:,:,:,1:split], x_sdf[2][:,:,:,1:split], 
        x_sdf[3][1:split], x_sdf[4][1:split]), y_puv[:,:,:,1:split]); 
        batchsize=16, shuffle=true, partial=false) |> device

    test_x = (x_sdf[1][:,:,:, split+1:end], x_sdf[2][:,:,:, split+1:end], 
        x_sdf[3][split+1:end], x_sdf[4][split+1:end])

    test_y = y_puv[:,:,:, split+1:end]

    return train_data, test_x, test_y
end

function load_data_small()
    @load "training_data/airfoil_training_data_small_aoa.jld2" x_sdf y_puv

    split = 64
    train_data = [((x_sdf[1][:,:,:,1:split], x_sdf[2][:,:,:,1:split], 
        x_sdf[3][1:split], x_sdf[4][1:split],
        x_sdf[5][1:split], x_sdf[6][1:split]), 
        1000y_puv[:,:,:,1:split])] |> device

    test_x = (x_sdf[1][:,:,:, split+1:end], x_sdf[2][:,:,:, split+1:end], 
        x_sdf[3][split+1:end], x_sdf[4][split+1:end],
        x_sdf[5][split+1:end], x_sdf[6][split+1:end]) |> device

    test_y = 1000y_puv[:,:,:, split+1:end] |> device

    return train_data, test_x, test_y
end

function plot_training_progress(x, y, n)
    cs =  heatmap(x[:,:,1,n]'; clim=(minimum(x[:,:,1,n]), maximum(x[:,:,1,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True Pressure")
    cui = heatmap(x[:,:,2,n]'; clim=(minimum(x[:,:,2,n]), maximum(x[:,:,2,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True U Velocity")
    cvi = heatmap(x[:,:,3,n]'; clim=(minimum(x[:,:,3,n]), maximum(x[:,:,3,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false,title="True V Velocity")
    cp =  heatmap(y[:,:,1,n]'; clim=(minimum(x[:,:,1,n]), maximum(x[:,:,1,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted Pressure")
    cu =  heatmap(y[:,:,2,n]'; clim=(minimum(x[:,:,2,n]), maximum(x[:,:,2,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted U Velocity")
    cv =  heatmap(y[:,:,3,n]'; clim=(minimum(x[:,:,3,n]), maximum(x[:,:,3,n])), size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted V Velocity")

    # clim=(minimum(y[:,:,1,n]), maximum(y[:,:,1,n])),
    # clim=(minimum(y[:,:,2,n]), maximum(y[:,:,2,n])),
    # clim=(minimum(y[:,:,3,n]), maximum(y[:,:,3,n])),
    # clim=(minimum(y[:,:,1,n]), maximum(y[:,:,1,n])),
    # clim=(minimum(y[:,:,2,n]), maximum(y[:,:,2,n])),
    # clim=(minimum(y[:,:,3,n]), maximum(y[:,:,3,n])),

    plot(cs, cp, cui, cu, cvi, cv; layout=(3,2))
end

function plot_training_progress_1(x, y, n)
    cs =  contour(x[:,:,1,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True Pressure")
    cui = contour(x[:,:,2,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True U Velocity")
    cvi = contour(x[:,:,3,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false,title="True V Velocity")
    cp =  contour(y[:,:,1,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted Pressure")
    cu =  contour(y[:,:,2,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted U Velocity")
    cv =  contour(y[:,:,3,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted V Velocity")

    plot(cs, cp, cui, cu, cvi, cv; layout=(3,2))
end

function plot_final(x, y, n)
    cs =  contour(x[:,:,1,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True Pressure")
    cui = contour(x[:,:,2,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="True U Velocity")
    cvi = contour(x[:,:,3,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false,title="True V Velocity")

    cp =  contour(y[:,:,1,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted Pressure")
    cu =  contour(y[:,:,2,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted U Velocity")
    cv =  contour(y[:,:,3,n]'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Predicted V Velocity")

    plot(cs, cp, cui, cu, cvi, cv; layout=(3,2))
end

function plot_diff(x, y, n)
    # dp =  contour(abs.(x[:,:,1,n]-y[:,:,1,n])'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Absolute P Difference")
    # du =  contour(abs.(x[:,:,2,n]-y[:,:,2,n])'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Absolute U Difference")
    # dv =  contour(abs.(x[:,:,3,n]-y[:,:,3,n])'; size=(1000,800), fill=(true, :jet), grid=false, xaxis=false, yaxis=false, title="Absolute V Difference")
end

function view_progress(x, y, u; n=1)
    ux = u(x)
    ux = ux .* x[2] |> cpu
    yy = y |> cpu

    plot_training_progress_1(yy, ux, n)
end

function view_final(x,y,u;n=1)
    ux = u(x)
    ux = ux .* x[2] |> cpu
    yy = y |> cpu

    yy = yy / std(yy) ./ x[2]
    ux = ux / std(ux) ./ x[2]

    plot_final(yy, ux, n)
end

function create_model_tangsali()
    act = leakyrelu
    depth = 128
    m_encode = Chain(
        Conv((5,5), 1 => depth, act; pad=(1,1)),
        MeanPool((4,4); pad=(1,1)),
    
        BatchNorm(depth),
        Conv((5,5), depth => depth, act; pad=(1,1)),
        MeanPool((4,4); pad=(1,1)),
    
        BatchNorm(depth),
        Conv((5,5), depth => depth, act; pad=(1,1)),
        MeanPool((4,4); pad=(1,1)),
    
        BatchNorm(depth),
        Conv((3,3), depth => depth, act; pad=(1,1)),
        Flux.flatten
    )
    
    m_decode = Chain(
        Dense(4depth + 4 => 4depth, act),
        a -> Flux.reshape(a, (2,2,depth,:)),
        Upsample(:bilinear, scale = (4,4)),
        
        BatchNorm(depth),
        Conv((3,3), depth => depth, act; pad=(1,1)),
        Upsample(:bilinear, scale = (4,4)),
    
        BatchNorm(depth),
        Conv((3,3), depth => depth, act; pad=(1,1)),
        Upsample(:bilinear, scale = (4,4)),
        BatchNorm(depth)
    )
    
    m_p = Conv((3,3), depth => 1, act; pad=(1,1))
    m_u = Conv((3,3), depth => 1, act; pad=(1,1))
    m_v = Conv((3,3), depth => 1, act; pad=(1,1))

    function model(x)
        y = [m_encode(x[1]) ; x[3]] |> m_decode
        p = y |> m_p
        u = y |> m_u
        v = y |> m_v
        return [p ;;; u ;;; v]
    end

    return model
end

function create_model_bhatnagar()
    act = swish
    depth = 128
    m_encode = Chain(
        Conv((5,5), 1 => depth, act; pad=(2,2)),
        MeanPool((4,4); pad=(1,1)),
    
        Conv((5,5), depth => depth, act; pad=(2,2)),
        MeanPool((4,4); pad=(1,1)),
    
        Conv((5,5), depth => depth, act; pad=(2,2)),
        MeanPool((4,4); pad=(1,1)),
    
        Conv((3,3), depth => depth, act; pad=(1,1)),
        Flux.flatten
    )
    
    m_decode = Chain(
        Dense(4depth + 2 => 4depth, act),
        a -> Flux.reshape(a, (2,2,depth,:)),
        Conv((3,3), depth => depth, act; pad=(1,1)),
        Upsample(:bilinear, scale = (4,4)),
        
        Conv((5,5), depth => depth, act; pad=(2,2)),
        Upsample(:bilinear, scale = (4,4)),
    
        Conv((5,5), depth => depth, act; pad=(2,2)),
        Upsample(:bilinear, scale = (4,4)),

        Conv((5,5), depth => 3, act; pad=(2,2)),
    )

    return function model(x)
        return [m_encode(x[1]) ; x[3]' ; x[4]'] |> m_decode
    end, Flux.params(m_encode, m_decode)
end

struct CustomModel
    encoder::Chain
    decoder::Chain
end

function CustomModel(;small=true)
    if small
        act = leakyrelu
        depth = 128
        m_encode = Chain(
            Conv((3,3), 1 => depth, act; pad=(1,1), bias=false),
            MeanPool((4,4); pad=(1,1)),
        
            Conv((3,3), depth => depth, act; pad=(1,1), bias=false),
            MeanPool((4,4); pad=(1,1)),
        
            Conv((3,3), depth => depth, act; pad=(1,1), bias=false),
            MeanPool((3,3); pad=(1,1)),
        
            Conv((3,3), depth => depth, act; pad=(1,1), bias=false),
            Flux.flatten
        ) |> device
        
        m_decode = Chain(
            Dense(4depth + 4 => 4depth, act, bias=false),
            a -> Flux.reshape(a, (2,2,depth,:)),
            Upsample(:bilinear, scale = (2,2)),
        
            Conv((3,3), depth => depth, act; pad=(1,1), bias=false),
            Upsample(:bilinear, scale = (4,4)),
            
            Conv((3,3), depth => depth, act; pad=(1,1), bias=false),
            Upsample(:bilinear, scale = (4,4)),
        
            Conv((3,3), depth => 3, identity; pad=(1,1), bias=false),
        ) |> device

        return CustomModel(m_encode, m_decode)
    else
        act = swish
        depth = 128
        m_encode = Chain(
            Conv((5,5), 1 => depth, act; pad=(2,2)),
            MeanPool((4,4); pad=(1,1)),
        
            Conv((5,5), depth => depth, act; pad=(2,2)),
            MeanPool((4,4); pad=(1,1)),
        
            Conv((5,5), depth => depth, act; pad=(2,2)),
            MeanPool((4,4); pad=(1,1)),
        
            Conv((3,3), depth => depth, act; pad=(1,1)),
            Flux.flatten
        ) |> device
        
        m_decode = Chain(
            Dense(4depth + 2 => 4depth, act),
            a -> Flux.reshape(a, (2,2,depth,:)),
            Conv((3,3), depth => depth, act; pad=(1,1)),
            Upsample(:bilinear, scale = (4,4)),
            
            Conv((5,5), depth => depth, act; pad=(2,2)),
            Upsample(:bilinear, scale = (4,4)),
        
            Conv((5,5), depth => depth, act; pad=(2,2)),
            Upsample(:bilinear, scale = (4,4)),

            Conv((5,5), depth => 3, act; pad=(2,2)),
        ) |> device

        return CustomModel(m_encode, m_decode)
    end
end

function (m::CustomModel)(x)
    return [m.encoder(x[1]) ; x[3]' ; x[4]' ; x[5]' ; x[6]'] |> m.decoder
end

Flux.@functor CustomModel |> device

function loss(x, y)
    mx = model(x)
    mask = x[2]
    l::Float32 = 0.0
    ϵ = 1e-5

    for i in axes(y,4)
        l += norm((mx[:,:,1:1,i] .* mask[:,:,1:1,i] - y[:,:,1:1,i]) / (std(y[:,:,1,i]) + ϵ)) + 
             norm((mx[:,:,2:2,i] .* mask[:,:,1:1,i] - y[:,:,2:2,i]) / (std(y[:,:,2,i]) + ϵ)) + 
             norm((mx[:,:,3:3,i] .* mask[:,:,1:1,i] - y[:,:,3:3,i]) / (std(y[:,:,3,i]) + ϵ))
    end

    d::Float32 = 0.0
    for i in 1:(size(y,4) - 1)
        d += ((std(y[:,:,1,i]) + std(y[:,:,1,i+1])) / 
             (norm(mx[:,:,1:1,i] .* mask[:,:,1:1,i] - 
             mx[:,:,1:1,i+1] .* mask[:,:,1:1,i+1]) + ϵ)) + 

             ((std(y[:,:,2,i]) + std(y[:,:,2,i+1])) / 
             (norm(mx[:,:,2:2,i] .* mask[:,:,1:1,i] - 
             mx[:,:,2:2,i+1] .* mask[:,:,1:1,i+1]) + ϵ)) + 

             ((std(y[:,:,3,i]) + std(y[:,:,3,i+1])) / 
             (norm(mx[:,:,3:3,i] .* mask[:,:,1:1,i] - 
             mx[:,:,3:3,i+1] .* mask[:,:,1:1,i+1]) + ϵ))
    end
    return (l + d) / prod(size(y[:,:,:,1]))
end

function loss_test(x, y)
    mx = x
    mask = x[2]
    l::Float32 = 0.0
    ϵ = 1e-5

    for i in size(y,4)
        l += norm((mx[:,:,1:1,i] - y[:,:,1:1,i]) / (std(y[:,:,1,i]) + ϵ)) + 
             norm((mx[:,:,2:2,i] - y[:,:,2:2,i]) / (std(y[:,:,2,i]) + ϵ)) + 
             norm((mx[:,:,3:3,i] - y[:,:,3:3,i]) / (std(y[:,:,3,i]) + ϵ))
    end
    return l / prod(size(y[:,:,:,1]))
end

function run()
    train_data, test_x, test_y = load_data_small()
    model = CustomModel(;small=true)
    opt = ADAM()

    test_losses = zeros(0)
    for epoch in 1:100
        train!(loss, Flux.params(model), train_data, opt)

        test_loss = loss(test_x, test_y)
        test_losses = [test_losses ; test_loss]

        @show i test_loss
    end
    return model
end

train_data, test_x, test_y = load_data_small();
model = CustomModel(;small=true);
opt = ADAM();

test_losses = zeros(0);
train_losses = zeros(0);
best_model = deepcopy(model)
@time for epoch in 1:150
    train!(loss, Flux.params(model), train_data, opt)

    global train_losses = [train_losses ; loss(train_data[1][1], train_data[1][2])]

    test_loss = loss(test_x, test_y)
    global test_losses = [test_losses ; test_loss]

    if test_loss < minimum(test_losses)
        global best_model = deepcopy(model)
    end
    @show epoch test_loss (test_y_norm - norm(model(test_x)))/test_y_norm
end

norm(model(test_x))
norm(test_y)

plot(test_losses; yaxis = :log)
view_progress(test_x, test_y, best_model; n=8)


@save "training_data/working_model_V6.jld2" model test_losses
test_losses[end]

loss_test(zeros(64,64,3,12) |> cpu, test_y |> cpu)

test_y |> size
train_data, test_x, test_y = load_data_small();
model = CustomModel(;small=true)
opt = ADAM();

# for d in train_data
#     @show model(d[1]) |> size
#     @show loss(d[1], d[2])
# end

# train!(loss, Flux.params(model), train_data, opt)

@load "training_data/working_model_V3.jld2" model test_losses

vf = view_final(test_x, test_y, model; n=2)
png(vf,"Paper/figures/model_v3.png")