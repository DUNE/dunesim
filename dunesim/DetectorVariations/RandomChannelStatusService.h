// Chris Backhouse - c.backhouse@ucl.ac.uk Dec 2017

#ifndef DETVAR_RANDOMCHANNELSTATUS_SERVICE
#define DETVAR_RANDOMCHANNELSTATUS_SERVICE

#include "larevt/CalibrationDBI/Interface/ChannelStatusService.h"

namespace detvar
{
  class RandomChannelStatusProvider: public lariov::ChannelStatusProvider
  {
  public:
    bool IsBad(raw::ChannelID_t chan) const override
    {
      return fBadChans.count(chan);
    }

    bool IsPresent(raw::ChannelID_t) const override {return true;}
    bool IsNoisy(raw::ChannelID_t) const override {return false;}

    std::set<raw::ChannelID_t> GoodChannels() const override
    {
      return fGoodChans;
    }

    std::set<raw::ChannelID_t> BadChannels() const override
    {
      return fBadChans;
    }

    std::set<raw::ChannelID_t> NoisyChannels() const override {return {};}

  protected:
    friend class RandomChannelStatusService;
    RandomChannelStatusProvider(const fhicl::ParameterSet& pset);

    std::set<raw::ChannelID_t> fBadChans, fGoodChans;
  };


  class RandomChannelStatusService: public lariov::ChannelStatusService
  {
  public:
    RandomChannelStatusService(const fhicl::ParameterSet& pset)
      : fProvider(pset)
    {
    }

  protected:
    const lariov::ChannelStatusProvider* DoGetProviderPtr() const override
    {
      return &fProvider;
    }

    const lariov::ChannelStatusProvider& DoGetProvider() const override
    {
      return fProvider;
    }

    RandomChannelStatusProvider fProvider;
  };

} // namespace

DECLARE_ART_SERVICE_INTERFACE_IMPL(detvar::RandomChannelStatusService, lariov::ChannelStatusService, LEGACY)

#endif
