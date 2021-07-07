from matplotlib import pyplot as plt



def get_figure_SE1E2IR_day(particle, params, n_clade, n_class):
    (fig, axes) = plt.subplots(3, 2, figsize=[10, 8])
    state_var = particle.statevar
    state_var = particle.statevar[particle.statevar[:, 0] > 0, :]      # if saved npz is not trimmed

    axes[0, 0].plot(state_var[:, 0], state_var[:, 1])
    for i in range(n_clade):
        print (i)
        axes[0, 1].plot(state_var[:, 0], state_var[:, i * n_class + 2] + state_var[:, i * n_class + 5], label="clade #" + str(i))
        axes[1, 0].plot(state_var[:, 0], state_var[:, i * n_class + 3] + state_var[:, i * n_class + 6], label="clade #" + str(i))
        axes[1, 1].plot(state_var[:, 0], state_var[:, i * n_class + 4] + state_var[:, i * n_class + 7], label="clade #" + str(i))
        axes[2, 0].plot(state_var[:, 0], state_var[:, i * n_class + 8], label="clade #" + str(i))

    axes[0, 0].set_ylabel('S')
    axes[0, 1].set_ylabel('E1')
    axes[1, 0].set_ylabel('E2')
    axes[1, 1].set_ylabel('I')
    axes[2, 0].set_ylabel('cumI')

    axes[0, 0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[0, 1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1, 0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1, 1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[2, 0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0, 0].set_ylim([0, params['N']])
    axes[2, 0].set_ylim([0, params['N']])

    axes[2, 0].legend();
    fig.tight_layout();

    return fig





def get_figure_SE1E2IR_day_oneclade(particle, params, str):
    (fig, axes) = plt.subplots(1, 2, figsize=[10, 4])
    state_var = particle.statevar
    #state_var = particle.statevar[particle.statevar[:, 0] > 0, :]      # if saved npz is not trimmed

    axes[0].plot(state_var[:, 0], state_var[:, 1], label = "S class")

    axes[1].plot(state_var[:, 0], state_var[:, 2] , label="E1 class" )
    axes[1].plot(state_var[:, 0], state_var[:, 3] + state_var[:, 5], label="E2 class" )
    axes[1].plot(state_var[:, 0], state_var[:, 4] + state_var[:, 6], label="I class" )
    axes[0].plot(state_var[:, 0], state_var[:, 7], label="cumulative I class")

    axes[0].set_ylabel('number of individuals')
    axes[1].set_ylabel('number of individuals')

    axes[0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0].set_ylim([0, params['N']])
    #axes[0, 1].set_ylim([0, params['N']])
    plt.title(str)
    axes[0].legend(); axes[1].legend();
    fig.tight_layout();

    return fig




def get_figure_SE1E2IR_day_oneclade_type2(particle, params, str):
    (fig, axes) = plt.subplots(1, 2, figsize=[10, 4])
    state_var = particle.statevar
    #state_var = particle.statevar[particle.statevar[:, 0] > 0, :]      # if saved npz is not trimmed

    axes[0].plot(state_var[:, 0], state_var[:, 1], label = "S class")

    axes[1].plot(state_var[:, 0], state_var[:, 2] , label="E1 class" )
    axes[1].plot(state_var[:, 0], state_var[:, 3] + state_var[:, 5]+ state_var[:, 4] + state_var[:, 6], label="E2 and I class" )
    axes[0].plot(state_var[:, 0], state_var[:, 7], label="cumulative I class")

    axes[0].set_ylabel('number of individuals')
    axes[1].set_ylabel('number of individuals')

    axes[0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0].set_ylim([0, params['N']])
    #axes[0, 1].set_ylim([0, params['N']])
    plt.title(str)
    axes[0].legend(); axes[1].legend();
    fig.tight_layout();

    return fig




def get_figure_SEIR_hetero_day_oneclade(particle, params, n_clade, n_class, str1):
    (fig, axes) = plt.subplots(1, 2, figsize=[10, 4])
    state_var = particle.statevar
    #state_var = particle.statevar[particle.statevar[:, 0] > 0, :]      # if saved npz is not trimmed

    axes[0].plot(state_var[:, 0], state_var[:, 1], label = "S class")
    axes[1].plot(state_var[:, 0], state_var[:, 2], label="E class" )
    axes[1].plot(state_var[:, 0], state_var[:, 3] + state_var[:, 4], label="I class" )
    axes[0].plot(state_var[:, 0], state_var[:, 5], label="cumulative I class")

    axes[0].set_ylabel('number of individuals')
    axes[1].set_ylabel('number of individuals')

    axes[0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0].set_ylim([0, params['N']])
    #axes[0, 1].set_ylim([0, params['N']])
    axes[0].legend(); axes[1].legend();

    plt.title(str1)


    fig.tight_layout();

    return fig




def get_figure_SEIR_day_oneclade(particle, params, str):
    (fig, axes) = plt.subplots(1, 2, figsize=[10, 4])
    state_var = particle.statevar
    #state_var = particle.statevar[particle.statevar[:, 0] > 0, :]      # if saved npz is not trimmed

    axes[0].plot(state_var[:, 0], state_var[:, 1], label = "S class")
    axes[1].plot(state_var[:, 0], state_var[:, 2], label="E class" )
    axes[1].plot(state_var[:, 0], state_var[:, 3], label="I class" )
    axes[0].plot(state_var[:, 0], state_var[:, 4], label="cumulative I class")

    axes[0].set_ylabel('number of individuals')
    axes[1].set_ylabel('number of individuals')

    axes[0].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])
    axes[1].set_xlim([min(state_var[:, 0]), max(state_var[:, 0])])

    axes[0].set_ylim([0, params['N']])
    #axes[0, 1].set_ylim([0, params['N']])

    axes[0].legend(); axes[1].legend();
    plt.title(str)
    fig.tight_layout();

    return fig
