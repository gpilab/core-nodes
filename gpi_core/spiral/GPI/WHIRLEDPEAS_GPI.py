#         NO CLINICAL USE.  THE SOFTWARE IS NOT INTENDED FOR COMMERCIAL
# PURPOSES AND SHOULD BE USED ONLY FOR NON-COMMERCIAL RESEARCH PURPOSES.  THE
# SOFTWARE MAY NOT IN ANY EVENT BE USED FOR ANY CLINICAL OR DIAGNOSTIC
# PURPOSES.  YOU ACKNOWLEDGE AND AGREE THAT THE SOFTWARE IS NOT INTENDED FOR
# USE IN ANY HIGH RISK OR STRICT LIABILITY ACTIVITY, INCLUDING BUT NOT LIMITED
# TO LIFE SUPPORT OR EMERGENCY MEDICAL OPERATIONS OR USES.  LICENSOR MAKES NO
# WARRANTY AND HAS NOR LIABILITY ARISING FROM ANY USE OF THE SOFTWARE IN ANY
# HIGH RISK OR STRICT LIABILITY ACTIVITIES.
#

# Author: Guru Krishnamoorthy
# Date: 2024Oct

import gpi
import numpy as np

class ExternalNode(gpi.NodeAPI):

    """Module to generate WHIRL (involute of a circle) trajectory and
    associated gradient waveforms using WHIRLED PEAS algorithm.  Also 
    generates associated SDC and time map for image construction.

    OUTPUTS:
    crds_out - output coordinates: the last dimension is 2 (kx/ky).
    grd_out - gradient waveforms used to produce crds_out.
    sdc_out - sampling density weights corresponding to crd array
    time_out - time map of k-space, for use in generating blurring kernels

    WIDGETS: self-explanatory
    Info Box gives the time for 4 WHIRL segments and respective constraints
    ** if constraint is less than what is requested, this was necessary to prevent a
       negative timing - it didn't slow anything down, it's just the max that was needed **
    """

    def execType(self):
        return gpi.GPI_PROCESS

    def initUI(self):
        
        # Widgets
        self.addWidget('PushButton', 'Compute', toggle=True)
        self.addWidget('TextBox', 'Info:')
        self.addWidget('ComboBox', 'Single Tau/Dual Tau', items=['Single', 'Dual'], visible=False)
        self.addWidget('ComboBox', 'Spiral Direction', items=['Spiral Out', 'Spiral In', 'Spiral In-Out', 'Spiral Out-In'])
        self.addWidget('DoubleSpinBox', 'FOV (cm)', val=24.0, min=0.1)
        self.addWidget('DoubleSpinBox', 'Res (mm)', val=0.8, min=0.1, singlestep=0.1)
        self.addWidget('DoubleSpinBox', 'Req. Tau for spiral-In (ms)', val=10.0, min=1, singlestep=0.1, visible=False)
        self.addWidget('DoubleSpinBox', 'Req. Tau for spiral-Out (ms)', val=10.0, min=1, singlestep=0.1, decimals=5)
        self.addWidget('DoubleSpinBox', 'Variable Density', val=0, min=0, max=2.0)
        self.addWidget('DoubleSpinBox', 'MaxSlw (mT/m/ms)', val=150.0, min=0.01)
        self.addWidget('DoubleSpinBox', 'MaxGrd (mT/m)', val=40.0, min=0.01)
        self.addWidget('DoubleSpinBox', 'Max G Freq (kHz)', val=1.0, min=0.1)
        self.addWidget('SpinBox', 'AD dwell time (us)', val=2, min=1)
        self.addWidget('SpinBox', 'Grad dwell time (us)', val=4, min=1)
        self.addWidget('DoubleSpinBox', 'Gam (kHz/mT)', val=42.577, min=0.01, visible=False)
        self.addWidget('SpinBox', 'Number fo Echoes', val=1, min=1, max=100)
        self.addWidget('ComboBox', 'Arm Order', items=['Linear', 'Skip', 'TwoWay', 'Mixed', 'Golden'])
        

        # Out Ports
        self.addOutPort('grd', 'NPYarray')
        self.addOutPort('crds', 'NPYarray')
        self.addOutPort('sdc', 'NPYarray')
        self.addOutPort('time_map', 'NPYarray')
        
        
        self.params = None


    def validate(self):
        
        spiral_direction = self.getVal('Spiral Direction')

        if spiral_direction == 'Spiral In-Out':
            self.setAttr('Single Tau/Dual Tau', visible=True)
            self.setAttr('Req. Tau for spiral-Out (ms)', visible=True)
            self.setAttr('Req. Tau for spiral-In (ms)', visible=self.getVal('Single Tau/Dual Tau') == 'Dual')
        elif spiral_direction == 'Spiral Out-In':
            self.setAttr('Single Tau/Dual Tau', visible=False)
            self.setAttr('Req. Tau for spiral-Out (ms)', visible=True)
            self.setAttr('Req. Tau for spiral-In (ms)', visible=False)
        elif spiral_direction == 'Spiral In':
            self.setAttr('Single Tau/Dual Tau', visible=False)
            self.setAttr('Req. Tau for spiral-Out (ms)', visible=False)
            self.setAttr('Req. Tau for spiral-In (ms)', visible=True)
        else:
            self.setAttr('Single Tau/Dual Tau', visible=False)
            self.setAttr('Req. Tau for spiral-In (ms)', visible=False)
            
        if self.getVal('Single Tau/Dual Tau') == 'Dual':
            self.setAttr('Number fo Echoes', val=1 , visible=False)
            self.setAttr('Spiral Direction', visible=False, val='Spiral In-Out')
        else:
            self.setAttr('Number fo Echoes', visible=True)
            self.setAttr('Spiral Direction', visible=True)

    def compute(self):

        import numpy as np
        
        ##########################
        # Compute waveforms
        ##########################
        if self.getVal('Compute'):
                
            params = {
                'fov': 0.01 * self.getVal('FOV (cm)'),
                'res': 0.001 * self.getVal('Res (mm)'),
                'req_tau_in': self.getVal('Req. Tau for spiral-In (ms)'),
                'req_tau_out': self.getVal('Req. Tau for spiral-Out (ms)'),
                'm_slew': self.getVal('MaxSlw (mT/m/ms)'),
                'm_grad': self.getVal('MaxGrd (mT/m)'),
                'm_omega': 2.*np.pi*self.getVal('Max G Freq (kHz)'),
                'dwell': 0.001 * self.getVal('AD dwell time (us)'),
                'grast': 0.001 * self.getVal('Grad dwell time (us)'),
                'gamma': self.getVal('Gam (kHz/mT)'),
                'nechoes': self.getVal('Number fo Echoes'),
                'spiral_direction': self.getVal('Spiral Direction'),
                'add_outer_ring': True, #self.getVal('Spiral Direction') == 'Spiral In-Out' and self.getVal('Number fo Echoes') > 1,
                'alpha_vd': self.getVal('Variable Density'),
                'nav_pts': 0,  # Assuming nav_pts is not used in single tau
                'arm_order': self.getVal('Arm Order')
            }

            if self.getVal('Single Tau/Dual Tau') == 'Dual' and self.getVal('Spiral Direction') == 'Spiral In-Out':
                g_out, crds, sdc, params_in, params_out, time_map = wpgen_dual_tau(params)
                params.update({'nyq_arms': min(params_in['nyq_arms'], params_out['nyq_arms'])})
            else:
                params.update({
                    'req_tau': params['req_tau_in'] if self.getVal('Spiral Direction') == 'Spiral In' else params['req_tau_out']
                })
                g_out, crds, sdc, params_out, time_map = wpgen_single_tau(params)
                params.update({'nyq_arms': params_out['nyq_arms']})
                
            [crds_full, sdc_full, grad_full] = gen_rotated_traj(crds, sdc, g_out,  params)    
            
            self.setData('crds', np.squeeze(crds_full))
            self.setData('sdc', np.squeeze(sdc_full))
            self.setData('grd', np.squeeze(grad_full))
            self.setData('time_map', time_map)
            
            info = f"********** Timing **********\n"
            info += f"Arc: {params_out['t_arc']:.2f} ms | Freq: {params_out['t_omega']:.2f} ms @ {params_out['omega_max']/(2.*np.pi):.2f} kHz\n"
            info += f"Slew: {params_out['t_slew']:.2f} ms @ {params_out['slew_max']:.2f} mT/m/ms | Grad: {params_out['t_grad']:.2f} ms @ {params_out['grad_max']:.2f} mT/m\n"

            if params['add_outer_ring']:
                info += f"Out Ring: {params_out['t_ring']:.2f} ms\n"

            info += f"Readout: {params_out['tau_total']:.2f} ms | Ramp: {params_out['t_ramp']:.2f} ms\n"
            info += f"Time Echo 1: {params_out['time_echo1']:.2f} ms"
            if params['nechoes'] > 1:
                info += f" | ΔTE: {params_out['delta_te']:.2f} ms"

            info += f"\n\n********** Sampling **********\n"
            info += f"Nyquist Arms: {params_out['nyq_arms']} | SNR Factor: {params_out['snr_factor']:.2f}\n"
            info += f"Grid: {params_out['grd_mtx']} x {params_out['grd_mtx']}\n"

            if 'params_in' in locals():
                info += f"\n********** Timing (In) **********\n"
                info += f"Arc: {params_in['t_arc']:.2f} ms | Freq: {params_in['t_omega']:.2f} ms @ {params_in['omega_max']/(2.*np.pi):.2f} kHz\n"
                info += f"Slew: {params_in['t_slew']:.2f} ms @ {params_in['slew_max']:.2f} mT/m/ms | Grad: {params_in['t_grad']:.2f} ms @ {params_in['grad_max']:.2f} mT/m\n"

                if params['add_outer_ring']:
                    info += f"Out Ring: {params_in['t_ring']:.2f} ms\n"

                info += f"Readout: {params_in['tau_total']:.2f} ms | Ramp: {params_in['t_ramp']:.2f} ms\n"
                info += f"Time Echo 1: {params_in['time_echo1']:.2f} ms"
                if params['nechoes'] > 1:
                    info += f" | ΔTE: {params_in['delta_te']:.2f} ms"

                info += f"\n\n********** Sampling (In) **********\n"
                info += f"Nyquist Arms: {params_in['nyq_arms']} | SNR Factor: {params_in['snr_factor']:.2f}\n"
                info += f"Grid: {params_in['grd_mtx']} x {params_in['grd_mtx']}\n"

            self.setAttr('Info:', val=info)
        
    
        return(0)


#%% # WPGEN dual Tau
def wpgen_dual_tau(params):
    
    import gpi_core.spiral.whirledpeas as whirledpeas
    
    fov = params['fov']
    res = params['res']
    req_tau_in = params['req_tau_in']
    req_tau_out = params['req_tau_out']
    m_slew = params['m_slew']
    m_grad = params['m_grad']
    m_omega = params['m_omega']
    grast = params['grast']
    dwell = params['dwell']
    gamma = params['gamma']
    add_outer_ring = params['add_outer_ring']
    nav_pts = int(params['nav_pts'] / int(dwell * 1000))
    
    wpgen_in = whirledpeas.WPGen(fov, res, req_tau_in, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring)
    
    wpgen_out = whirledpeas.WPGen(fov, res, req_tau_out, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring)
    
    wpcompose = whirledpeas.WPCompose(wpgen_in, wpgen_out, nav_pts)
    
    g_out = wpcompose.composeGradients()
    [ksp, sdc] = wpcompose.composeKSPnSDC()
    
    [tmap_in, tmap_out] = wpcompose.composeTimeMap()
    time_map = None
    if np.any(~np.isnan(tmap_in)) and np.any(~np.isnan(tmap_out)):
        time_map = np.stack((tmap_in, tmap_out), axis=0)
    elif np.any(~np.isnan(tmap_in)):
        time_map = tmap_in
    elif np.any(~np.isnan(tmap_out)):
        time_map = tmap_out
    
    def extract_params(wpgen, wpcompose):
        return {
            'time_echo1': wpcompose.getTimeEcho1(),
            'delta_te': wpcompose.getDeltaTE(),
            'grad_max': wpgen.GetGradMax(),
            'nyq_arms': wpgen.GetNyqNumberArms(),
            't_arc': wpgen.GetArcDur(),
            't_omega': wpgen.GetFreqDur(),
            't_slew': wpgen.GetSlewDur(),
            't_grad': wpgen.GetGradDur(),
            'tau_total': wpgen.GetActualTau(),
            't_ramp': wpgen.GetRampDur(),
            't_ring': wpgen.GetOutRingDur(),
            'omega_max': wpgen.GetOmegaMax(),
            'slew_max': wpgen.GetSlewMax(),
            'snr_factor': wpgen.GetSnrFactor(),
            'ksp_pts': wpgen.GetKspacePts(),
            'grd_mtx': wpgen.GetGridMtxSize()
        }

    params_in = extract_params(wpgen_in, wpcompose)
    params_out = extract_params(wpgen_out, wpcompose)
    
    return g_out, ksp, sdc, params_in, params_out, time_map

#%% # WPGEN single Tau
def wpgen_single_tau(params):
    
    import gpi_core.spiral.whirledpeas as whirledpeas
    
    fov = params['fov']
    res = params['res']
    req_tau = params['req_tau']
    m_slew = params['m_slew']
    m_grad = params['m_grad']
    m_omega = params['m_omega']
    grast = params['grast']
    dwell = params['dwell']
    gamma = params['gamma']
    add_outer_ring = params['add_outer_ring']
    nechoes = params['nechoes']
    alpha_vd = params['alpha_vd']
    spiral_direct = {
                    'Spiral Out': whirledpeas.SpiralDirect.SPIRAL_OUT,
                    'Spiral In': whirledpeas.SpiralDirect.SPIRAL_IN,
                    'Spiral In-Out': whirledpeas.SpiralDirect.SPIRAL_INOUT,
                    'Spiral Out-In': whirledpeas.SpiralDirect.SPIRAL_OUTIN
                }.get(params['spiral_direction'], whirledpeas.SpiralDirect.SPIRAL_OUT)
    
    wpgen = whirledpeas.WPGen(fov, res, req_tau, m_slew, m_grad, m_omega, grast, dwell, gamma,
                                    add_outer_ring, alpha_vd)
    
    wpcompose = whirledpeas.WPCompose(wpgen, nechoes, spiral_direct)
    
    g_out = wpcompose.composeGradients()
    [ksp, sdc] = wpcompose.composeKSPnSDC()
    [tmap_in, tmap_out] = wpcompose.composeTimeMap()
    time_map = None
    if np.any(~np.isnan(tmap_in)) and np.any(~np.isnan(tmap_out)):
        time_map = np.stack((tmap_in, tmap_out), axis=0)
    elif np.any(~np.isnan(tmap_in)):
        time_map = tmap_in
    elif np.any(~np.isnan(tmap_out)):
        time_map = tmap_out
    
    out_params = {
        'time_echo1': wpcompose.getTimeEcho1(),
        'delta_te': wpcompose.getDeltaTE(),
        't_arc': wpgen.GetArcDur(),
        't_omega': wpgen.GetFreqDur(),
        't_slew': wpgen.GetSlewDur(),
        't_grad': wpgen.GetGradDur(),
        'tau_total': wpgen.GetActualTau(),
        't_ramp': wpgen.GetRampDur(),
        't_ring': wpgen.GetOutRingDur(),
        'omega_max': wpgen.GetOmegaMax(),
        'slew_max': wpgen.GetSlewMax(),
        'grad_max': wpcompose.getGradMax(),
        'snr_factor': wpgen.GetSnrFactor(),
        'ksp_pts': wpgen.GetKspacePts(),
        'nyq_arms': wpcompose.getNyqArms(),
        'grd_mtx': wpgen.GetGridMtxSize()
    }
    
    return g_out, ksp, sdc, out_params, time_map

def gen_rotated_traj(ksp, sdc, grad, params):
    
    import gpi_core.spiral.whirledpeas as whirledpeas
    
    arm_ord = {
                'Linear': whirledpeas.OrderType.LINEAR,
                'Skip': whirledpeas.OrderType.SKIP,
                'TwoWay': whirledpeas.OrderType.TWO_WAY,
                'Mixed': whirledpeas.OrderType.MIXED,
                'Golden': whirledpeas.OrderType.GOLDEN
            }.get(params['arm_order'], whirledpeas.OrderType.LINEAR)
    
    angular_coverage = np.pi if params['spiral_direction'] == 'Spiral In-Out' else 2*np.pi
    sporder = whirledpeas.SpiralArmOrder(params['nyq_arms'],  1, angular_coverage,arm_ord,  False, False)
    arm_angles = sporder.composeAngles()

    crd_full = rotate_trajectory(ksp, arm_angles, params['nyq_arms'], params['nechoes'])
    sdc_full = np.broadcast_to(sdc, (params['nechoes'], params['nyq_arms'], sdc.shape[-1]))
    sdc_full = np.squeeze(sdc_full)
    
    grad_full = rotate_trajectory(grad, arm_angles, params['nyq_arms'], 1)
    
    return crd_full, sdc_full, grad_full

#%% # Rotate k-space trajectory
def rotate_trajectory(arm_to_rotate, arm_angles, total_arms, total_echo):
    sample_pts = arm_to_rotate.shape[-2]
    coords_full = np.zeros([total_echo, total_arms, sample_pts, 2], dtype=np.double)
    if arm_to_rotate.ndim == 2:
        arm_to_rotate = np.expand_dims(arm_to_rotate, axis=-3)
        
    for arm in range(total_arms):
        for echo in range(total_echo):
            curr_angle = arm_angles[0, arm]
            sb = np.sin(curr_angle)
            cb = np.cos(curr_angle)
            coords_full[echo, arm, :, 0] = (cb * arm_to_rotate[echo, :, 0]) - (sb * arm_to_rotate[echo, :, 1])
            coords_full[echo, arm, :, 1] = (cb * arm_to_rotate[echo, :, 1]) + (sb * arm_to_rotate[echo, :, 0])
    
    coords_full = np.squeeze(coords_full)
    return coords_full
